% =========================================================================
% DoA_KW_UCA_SelfCal_MultiDir
%
% Self-calibration iterativa por COMPENSACAO de acoplamento usando MULTIPLAS
% DIRECOES NAO-SIMULTANEAS, com foco em MELHORAR e DIAGNOSTICAR a estimacao
% da matriz C.
%
%   Modelo (P medidas, nao-simultaneas):
%       b_p = alpha_p * C * a(phi_p) + ruido,   p = 1..P    (C circulante)
%
%   Laco alternante (init = KW):
%     1) para cada p: phi_p = DoA( inv(C) * X_p )
%     2) C <- estimate_C_md( {b_p}, {a(phi_p)} )   (LS conjunto, opc. ridge)
%     3) repete
%
%   MELHORIAS desta versao (flags no topo):
%     - drive_C_with_KW : a DIRECAO que alimenta a estimacao de C vem SEMPRE
%                         do KW (robusto a acoplamento), evitando a propagacao
%                         de erro de DoA / ciclos-limite dos metodos espectrais.
%                         A DoA REPORTADA continua sendo a do metodo avaliado.
%     - ridge_lambda    : Tikhonov no LS dos coeficientes (doma o mau
%                         condicionamento em aberturas pequenas; 0 = desliga).
%     - stratified_dirs : sorteio das P direcoes por faixas (cobertura angular
%                         garantida), melhorando a identificabilidade.
%     - safe_inv        : inversao robusta (evita NaN/Inf em C mal-condicionada).
%
%   METRICAS de C (alem do Frobenius bruto, que e' enganoso):
%     - Frobenius invariante a escala : min_s ||C_hat - s C_true|| / ||C_true||
%     - Residuo de compensacao efetiva: ||norm(inv(C_hat) C_true) - I|| / sqrt(M)
%       (mede o que importa: quao bem inv(C_hat) DESFAZ o acoplamento)
%     - Erro por coeficiente |c_hat_k - c_k| e PLANO COMPLEXO (verdadeiro vs
%       nuvem de estimativas), pois ha so K = floor(M/2) coef. unicos.
% =========================================================================

clear; clc; close all;

%% ---- Parametros do array ----
M      = 8;            % nº de elementos do UCA
fc     = 500e6;        % Hz
c      = 3e8;
lambda = c/fc;
r      = 0.15*lambda;   % raio (0.1 lambda = abertura pequena, acoplamento forte)
theta_sig_deg = 90;    % elevacao (plano XY)

%% ---- Parametros do sinal (FSK-2 conhecido) ----
% NOTA: aumentar N (preambulo conhecido) reduz a variancia de b_hat e melhora C.
fs = 288000; N = 3000; Rs = 9600; sps = 30; alpha = 0.3; span = 8; fd = 4.8e3;

%% ---- Parametros do experimento ----
maxIter      = 6;
range_SNR_dB = -12:3:12; %[-6, 0, 6, 12];
methods      = {'KW','DAS','CAPON','MUSIC'};
nMethods     = numel(methods);
phi_grid_deg = -180:0.5:180;
grid_step    = phi_grid_deg(2) - phi_grid_deg(1);
beta_uca     = 2*pi*(0:M-1).'/M;
K            = floor(M/2);          % nº de coeficientes circulantes unicos

% --- Multi-direcao ---
P_default    = 5;
range_P      = [1 2 3 4 6 8];
snr_Psweep   = 6;
min_sep_deg  = 10;
inner_iter   = 30;
inner_tol    = 1e-10;

% --- MELHORIAS (flags) ---
drive_C_with_KW = false;     % direcoes que alimentam C vem do KW (robusto)
ridge_lambda    = 5e-3;     % Tikhonov no LS de c (0 = sem regularizacao)
stratified_dirs = true;     % cobertura angular garantida no sorteio das P direcoes

% --- Monte Carlo ---
n_runs   = 200;
rng(2026, 'twister');

if drive_C_with_KW, drv_tag = 'C dirigida por KW'; else, drv_tag = 'C dirigida pelo metodo'; end
outDir = fullfile(pwd, 'MultiDir SelfCal Graphs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---- Matriz de acoplamento VERDADEIRA (so' p/ medir erro) ----
C_true    = compute_Ctx_for_R(fc, M, r, 50);
normCtrue = norm(C_true, 'fro');
c_true_vec = C_true(1, 2:K+1).';     % coeficientes unicos verdadeiros (K x 1)
fprintf('C_true (raio=%.2f lambda, ||C_true||_F=%.4f, cond=%.1f)\n', ...
    r/lambda, normCtrue, cond(C_true));

%% ---- Dicionario de steering vectors ----
nGrid = numel(phi_grid_deg);
A_dict = zeros(M, nGrid);
for ig = 1:nGrid
    A_dict(:,ig) = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_grid_deg(ig));
end

%% =====================================================================
%   PARTE A/C: convergencia por SNR e resumo vs SNR  (P = P_default)
%% =====================================================================
nSNR             = numel(range_SNR_dB);
sumsq_doa_iter   = zeros(nMethods, maxIter, nSNR);
sum_frob_iter    = zeros(nMethods, maxIter, nSNR);
sum_compres_iter = zeros(nMethods, maxIter, nSNR);
sumsq_doa_noComp = zeros(nMethods, nSNR);
sumsq_doa_noMC   = zeros(nMethods, nSNR);
cnt_doa          = zeros(1, nSNR);
cnt_run          = zeros(1, nSNR);
c_runs           = zeros(nMethods, K, n_runs, nSNR);   % coef. estimados (final iter)

t0 = tic;
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fprintf('\n===== [A/C] SNR = %+d dB  (P=%d, %d rodadas, %s) =====\n', ...
        SNR_dB, P_default, n_runs, drv_tag);
    for run = 1:n_runs
        [B, Xc_cell, q_cell, phitrue] = gen_shots(P_default, grid_step, min_sep_deg, ...
            stratified_dirs, M, r, lambda, theta_sig_deg, SNR_dB, N, fs, Rs, sps, ...
            alpha, span, fd, C_true);

        % referencias (sem iterar)
        for im = 1:nMethods
            for p = 1:P_default
                phiMC = doa_estimate(q_cell{p}.X_ideal, methods{im}, q_cell{p}.q, r, lambda, beta_uca, A_dict, phi_grid_deg);
                phiNC = doa_estimate(Xc_cell{p},        methods{im}, q_cell{p}.q, r, lambda, beta_uca, A_dict, phi_grid_deg);
                sumsq_doa_noMC(im,iSNR)   = sumsq_doa_noMC(im,iSNR)   + wrapTo180(phiMC - phitrue(p))^2;
                sumsq_doa_noComp(im,iSNR) = sumsq_doa_noComp(im,iSNR) + wrapTo180(phiNC - phitrue(p))^2;
            end
        end

        % laco multi-dir por metodo
        for im = 1:nMethods
            dm = methods{im};  if drive_C_with_KW, dm = 'KW'; end
            [frob_it, compres_it, ssq_it, c_fin] = multidir_selfcal_one(methods{im}, dm, ...
                B, Xc_cell, q_cell, phitrue, M, r, lambda, theta_sig_deg, beta_uca, ...
                A_dict, phi_grid_deg, maxIter, C_true, normCtrue, inner_iter, inner_tol, ridge_lambda);
            sum_frob_iter(im,:,iSNR)    = sum_frob_iter(im,:,iSNR)    + frob_it(:).';
            sum_compres_iter(im,:,iSNR) = sum_compres_iter(im,:,iSNR) + compres_it(:).';
            sumsq_doa_iter(im,:,iSNR)   = sumsq_doa_iter(im,:,iSNR)   + ssq_it(:).';
            c_runs(im,:,run,iSNR)       = c_fin(:).';
        end
        cnt_doa(iSNR) = cnt_doa(iSNR) + P_default;
        cnt_run(iSNR) = cnt_run(iSNR) + 1;
    end
    fprintf('  decorrido %.0fs\n', toc(t0));
end

% --- Consolidacao A/C ---
rmse_doa_iter  = zeros(nMethods, maxIter, nSNR);
mean_frob_iter = zeros(nMethods, maxIter, nSNR);
mean_cres_iter = zeros(nMethods, maxIter, nSNR);
for iSNR = 1:nSNR
    rmse_doa_iter(:,:,iSNR)  = sqrt(sumsq_doa_iter(:,:,iSNR) / cnt_doa(iSNR));
    mean_frob_iter(:,:,iSNR) =      sum_frob_iter(:,:,iSNR)  / cnt_run(iSNR);
    mean_cres_iter(:,:,iSNR) =      sum_compres_iter(:,:,iSNR)/ cnt_run(iSNR);
end
rmse_doa_final = squeeze(rmse_doa_iter(:,maxIter,:));
frob_final     = squeeze(mean_frob_iter(:,maxIter,:));
cres_final     = squeeze(mean_cres_iter(:,maxIter,:));
rmse_doa_noComp= sqrt(sumsq_doa_noComp ./ cnt_doa);
rmse_doa_noMC  = sqrt(sumsq_doa_noMC   ./ cnt_doa);

% Metricas finais derivadas dos coeficientes (Frobenius invariante a escala
% e erro por coeficiente), medias sobre rodadas.
scaleinv_final = zeros(nMethods, nSNR);
coeff_err      = zeros(nMethods, K, nSNR);
for iSNR = 1:nSNR
    for im = 1:nMethods
        acc = 0;
        for run = 1:n_runs
            Chat = reconstruct_C(squeeze(c_runs(im,:,run,iSNR)).', M);
            acc  = acc + scaleinv_frob(Chat, C_true);
        end
        scaleinv_final(im,iSNR) = acc / n_runs;
        for k = 1:K
            coeff_err(im,k,iSNR) = mean(abs(squeeze(c_runs(im,k,:,iSNR)) - c_true_vec(k)));
        end
    end
end

flr = @(x) max(x, 1e-3);
colorsM   = lines(nMethods);
markers_m = {'o-','s-','d-','^-'};
% Indices de metodos a exibir nas figuras de METRICA DE C: se C e' dirigida
% por KW, todos os metodos compartilham a MESMA C -> mostra so uma curva.
if drive_C_with_KW, Cidx = 1; Clab = {'C (KW-driven)'}; else, Cidx = 1:nMethods; Clab = methods; end

% --- (A) convergencia por SNR: DoA, Frobenius e residuo de compensacao ---
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fig = figure('Name',sprintf('MultiDir SNR=%+d dB',SNR_dB),'Color','w','Position',[40 60 1500 470]);

    subplot(1,3,1); hold on; grid on;
    for im = 1:nMethods
        plot(1:maxIter, flr(squeeze(rmse_doa_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:),'LineWidth',1.6,'MarkerFaceColor',colorsM(im,:), ...
            'DisplayName',methods{im});
        yline(flr(rmse_doa_noComp(im,iSNR)), ':','Color',colorsM(im,:),'LineWidth',1.0,'HandleVisibility','off');
    end
    set(gca,'YScale','log'); xlabel('Iteracao'); ylabel('RMSE de DoA (graus)');
    title('DoA (pontilhado = sem comp.)'); legend('Location','best'); xlim([1 maxIter]);

    subplot(1,3,2); hold on; grid on;
    for jj = 1:numel(Cidx)
        im = Cidx(jj);
        plot(1:maxIter, flr(squeeze(mean_frob_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:),'LineWidth',1.6,'MarkerFaceColor',colorsM(im,:),'DisplayName',Clab{jj});
        plot(1:maxIter, flr(squeeze(mean_cres_iter(im,:,iSNR))), [markers_m{im}(1) '--'], ...
            'Color',colorsM(im,:),'LineWidth',1.2,'MarkerFaceColor','none', ...
            'DisplayName',[Clab{jj} ' (comp.res)']);
    end
    set(gca,'YScale','log'); xlabel('Iteracao'); ylabel('erro de C');
    title('Frobenius (solido) vs residuo de compensacao (tracejado)');
    legend('Location','best'); xlim([1 maxIter]);

    subplot(1,3,3); hold on; grid on;
    for k = 1:K
        re = real(c_true_vec(k)); imv = imag(c_true_vec(k));
        plot(re, imv, 'ko', 'MarkerSize',9, 'LineWidth',1.4, 'HandleVisibility','off');
        text(re, imv, sprintf(' c_%d',k), 'FontWeight','bold');
        re_e = squeeze(real(c_runs(Cidx(1),k,:,iSNR)));
        im_e = squeeze(imag(c_runs(Cidx(1),k,:,iSNR)));
        scatter(re_e, im_e, 12, 'filled', 'MarkerFaceAlpha',0.35, 'DisplayName',sprintf('c_%d est',k));
    end
    axis equal; grid on; xlabel('Re'); ylabel('Im');
    title(sprintf('Coeficientes (%s)', Clab{1}));

    sgtitle(sprintf('MultiDir  SNR=%+d dB,  P=%d,  raio=%.2f\\lambda,  %s,  ridge=%.0e', ...
        SNR_dB, P_default, r/lambda, drv_tag, ridge_lambda), 'FontWeight','bold');
    exportgraphics(fig, fullfile(outDir, sprintf('multidir_iter_SNR_%+03d.png',SNR_dB)),'Resolution',170);
    matlab2tikz(fullfile(outDir, sprintf('multidir_iter_SNR_%+03d.tex',SNR_dB)), 'width','\figurewidth','height','\figureheight');
end

% --- (C) resumo DoA vs SNR ---
fig = figure('Color','w','Position',[100 100 920 560]); hold on; grid on;
for im = 1:nMethods
    plot(range_SNR_dB, flr(rmse_doa_final(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',2.0,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',sprintf('%s multi-dir',methods{im}));
    plot(range_SNR_dB, flr(rmse_doa_noComp(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.2,'LineStyle',':','MarkerFaceColor','none','MarkerSize',8,'DisplayName',sprintf('%s sem comp.',methods{im}));
end
set(gca,'YScale','log'); xlabel('SNR (dB)'); ylabel('RMSE de DoA (graus)');
xticks(range_SNR_dB); xlim([min(range_SNR_dB)-1, max(range_SNR_dB)+1]);
title(sprintf('RMSE de DoA: multi-dir (solido) vs sem comp (pontilhado)\nP=%d, %d rodadas, raio=%.2f\\lambda, %s', ...
    P_default, n_runs, r/lambda, drv_tag));
legend('Location','eastoutside');
exportgraphics(fig, fullfile(outDir,'resumo_RMSE_doa_vs_snr.png'),'Resolution',180);
matlab2tikz(fullfile(outDir,'resumo_RMSE_doa_vs_snr.tex'), 'width','\figurewidth','height','\figureheight');

% --- (D) resumo METRICAS DE C vs SNR (Frobenius vs invariante a escala vs comp.res) ---
fig = figure('Color','w','Position',[60 80 1500 470]);
metsD = {frob_final, scaleinv_final, cres_final};
titD  = {'Frobenius bruto  ||C_{hat}-C_{true}||/||C_{true}||', ...
         'Frobenius invariante a escala', ...
         'Residuo de compensacao  ||inv(C_{hat})C_{true}-I||/\surd M'};
for sp = 1:3
    subplot(1,3,sp); hold on; grid on;
    for jj = 1:numel(Cidx)
        im = Cidx(jj);
        plot(range_SNR_dB, flr(metsD{sp}(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
            'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',Clab{jj});
    end
    set(gca,'YScale','log'); xlabel('SNR (dB)'); ylabel('erro de C'); xticks(range_SNR_dB);
    title(titD{sp}); legend('Location','best');
end
sgtitle(sprintf('Metricas de C vs SNR  (P=%d, %d rodadas, raio=%.2f\\lambda, %s, ridge=%.0e)', ...
    P_default, n_runs, r/lambda, drv_tag, ridge_lambda),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'resumo_metricas_C_vs_snr.png'),'Resolution',170);
matlab2tikz(fullfile(outDir,'resumo_metricas_C_vs_snr.tex'), 'width','\figurewidth','height','\figureheight');

% --- (E) PLANO COMPLEXO dos coeficientes, por SNR ---
ck_colors = [0.85 0.10 0.10; 0.10 0.60 0.15; 0.10 0.30 0.85; 0.60 0.15 0.65];
for iSNR = 1:nSNR
    if drive_C_with_KW, midx = 1; rows = 1; cols = 1; figpos=[200 200 640 560];
    else,               midx = 1:nMethods; rows = 2; cols = 2; figpos=[80 60 1100 820]; end
    fig = figure('Name',sprintf('Coef. complexos SNR=%+d',range_SNR_dB(iSNR)),'Color','w','Position',figpos);
    for jj = 1:numel(midx)
        im = midx(jj);
        subplot(rows,cols,jj); hold on; grid on; axis equal;
        for k = 1:K
            re_e = squeeze(real(c_runs(im,k,:,iSNR)));
            im_e = squeeze(imag(c_runs(im,k,:,iSNR)));
            scatter(re_e, im_e, 16, ck_colors(k,:), 'filled', 'MarkerFaceAlpha',0.30, ...
                'DisplayName',sprintf('c_%d (est)',k));
        end
        for k = 1:K
            plot(real(c_true_vec(k)), imag(c_true_vec(k)), 'p', 'MarkerSize',16, ...
                'MarkerFaceColor',ck_colors(k,:), 'MarkerEdgeColor','k', 'LineWidth',1.3, 'HandleVisibility','off');
            text(real(c_true_vec(k)), imag(c_true_vec(k)), sprintf(' c_%d',k), 'FontWeight','bold','FontSize',10);
        end
        xlabel('Re'); ylabel('Im');
        if drive_C_with_KW, ttl = 'C (KW-driven)'; else, ttl = methods{im}; end
        estr = strjoin(arrayfun(@(k) sprintf('%.2f',coeff_err(im,k,iSNR)), 1:K, 'uni',0), ', ');
        title(sprintf('%s   |\\Deltac_k|=[%s]', ttl, estr), 'FontSize',9);
        if jj == 1, legend('Location','bestoutside'); end
    end
    sgtitle(sprintf(['Coeficientes c_1..c_%d: verdadeiro (estrela preta) vs estimado (nuvem)\n' ...
        'SNR=%+d dB, %d rodadas, raio=%.2f\\lambda, %s'], ...
        K, range_SNR_dB(iSNR), n_runs, r/lambda, drv_tag), 'FontWeight','bold');
    exportgraphics(fig, fullfile(outDir, sprintf('coeficientes_complexo_SNR_%+03d.png',range_SNR_dB(iSNR))),'Resolution',170);
    matlab2tikz(fullfile(outDir, sprintf('coeficientes_complexo_SNR_%+03d.tex',range_SNR_dB(iSNR))), 'width','\figurewidth','height','\figureheight');
end

%% =====================================================================
%   PARTE B: varredura do numero de direcoes P  (SNR = snr_Psweep)
%% =====================================================================
fprintf('\n===== [B] Varredura de P em SNR=%+d dB =====\n', snr_Psweep);
nP = numel(range_P);
finalRMSE_vsP = zeros(nMethods, nP);
finalFrob_vsP = zeros(nMethods, nP);
finalCres_vsP = zeros(nMethods, nP);
for iP = 1:nP
    P = range_P(iP);
    sumsq_doa_P=zeros(nMethods,1); sum_frob_P=zeros(nMethods,1); sum_cres_P=zeros(nMethods,1);
    cnt_doa_P=0; cnt_run_P=0;
    for run = 1:n_runs
        [B, Xc_cell, q_cell, phitrue] = gen_shots(P, grid_step, min_sep_deg, ...
            stratified_dirs, M, r, lambda, theta_sig_deg, snr_Psweep, N, fs, Rs, sps, ...
            alpha, span, fd, C_true);
        for im = 1:nMethods
            dm = methods{im};  if drive_C_with_KW, dm = 'KW'; end
            [frob_it, compres_it, ssq_it, ~] = multidir_selfcal_one(methods{im}, dm, ...
                B, Xc_cell, q_cell, phitrue, M, r, lambda, theta_sig_deg, beta_uca, ...
                A_dict, phi_grid_deg, maxIter, C_true, normCtrue, inner_iter, inner_tol, ridge_lambda);
            sumsq_doa_P(im)=sumsq_doa_P(im)+ssq_it(end);
            sum_frob_P(im) =sum_frob_P(im) +frob_it(end);
            sum_cres_P(im) =sum_cres_P(im) +compres_it(end);
        end
        cnt_doa_P=cnt_doa_P+P; cnt_run_P=cnt_run_P+1;
    end
    finalRMSE_vsP(:,iP)=sqrt(sumsq_doa_P/cnt_doa_P);
    finalFrob_vsP(:,iP)=sum_frob_P/cnt_run_P;
    finalCres_vsP(:,iP)=sum_cres_P/cnt_run_P;
    fprintf('  P=%d feito (%.0fs)\n', P, toc(t0));
end

fig = figure('Color','w','Position',[60 80 1500 470]);
metsB = {finalRMSE_vsP, finalFrob_vsP, finalCres_vsP};
titB  = {'RMSE de DoA final (graus)', 'Frobenius bruto final', 'Residuo de compensacao final'};
idxB  = {1:nMethods, Cidx, Cidx};   % DoA por metodo; C-metricas conforme Cidx
labB  = {methods, Clab, Clab};
for sp = 1:3
    subplot(1,3,sp); hold on; grid on;
    ids = idxB{sp};
    for jj = 1:numel(ids)
        im = ids(jj);
        plot(range_P, flr(metsB{sp}(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
            'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',labB{sp}{jj});
    end
    set(gca,'YScale','log'); xlabel('nº de direcoes P'); ylabel(titB{sp}); xticks(range_P);
    title(titB{sp}); legend('Location','best');
end
sgtitle(sprintf('Ganho com mais direcoes  (SNR=%+d dB, %d rodadas, raio=%.2f\\lambda, %s)', ...
    snr_Psweep, n_runs, r/lambda, drv_tag),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'varredura_numero_direcoes_P.png'),'Resolution',170);
matlab2tikz(fullfile(outDir,'varredura_numero_direcoes_P.tex'), 'width','\figurewidth','height','\figureheight');

%% ---- Resumo numerico ----
fprintf('\n===== RESUMO (P=%d, %s, ridge=%.0e) =====\n', P_default, drv_tag, ridge_lambda);
for iSNR = 1:nSNR
    fprintf('SNR=%+3d dB:\n', range_SNR_dB(iSNR));
    for jj = 1:numel(Cidx)
        im = Cidx(jj);
        fprintf('   %-12s  Frob=%.3f  FrobEsc=%.3f  compRes=%.3f  |dc_k|=[%s]\n', ...
            Clab{jj}, frob_final(im,iSNR), scaleinv_final(im,iSNR), cres_final(im,iSNR), ...
            strjoin(arrayfun(@(k) sprintf('%.2f',coeff_err(im,k,iSNR)),1:K,'uni',0),', '));
    end
    for im = 1:nMethods
        fprintf('     DoA %-6s  multiDir=%7.3f  semComp=%7.3f\n', ...
            methods{im}, rmse_doa_final(im,iSNR), rmse_doa_noComp(im,iSNR));
    end
end
fprintf('\nFiguras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function [B, Xc_cell, q_cell, phitrue] = gen_shots(P, grid_step, min_sep, stratified, ...
    M, r, lambda, theta_deg, SNR_dB, N, fs, Rs, sps, alpha, span, fd, C_true)
% GEN_SHOTS  Sorteia P direcoes distintas e gera P shots NAO-SIMULTANEOS.
    phitrue = draw_dirs(P, grid_step, min_sep, stratified);
    B = zeros(M, P); Xc_cell = cell(1,P); q_cell = cell(1,P);
    for p = 1:P
        [~, q, ~, Xsig, ~, Xn] = utils.simulate_fsk_data_uca(M, r, lambda, ...
            phitrue(p), phitrue(p)+90, theta_deg, theta_deg, ...
            SNR_dB, -100, N, fs, Rs, sps, alpha, span, fd);
        q = q(:);
        Xc = C_true*Xsig + Xn;
        Xc_cell{p} = Xc;
        s.q = q; s.X_ideal = Xsig + Xn;
        q_cell{p} = s;
        B(:,p) = Xc*conj(q)/(q'*q);
    end
end

function ph = draw_dirs(P, grid_step, min_sep, stratified)
% DRAW_DIRS  P azimutes no grid; se stratified, um por faixa (cobertura ampla).
    ph = zeros(1, P);
    if stratified
        edges = linspace(-180, 180, P+1);
        for p = 1:P
            tries = 0;
            while true
                cand = edges(p) + (edges(p+1)-edges(p))*rand;
                cand = wrapTo180(round(cand/grid_step)*grid_step);
                if p == 1 || all(abs(wrapTo180(cand - ph(1:p-1))) >= min_sep)
                    ph(p) = cand; break;
                end
                tries = tries + 1;
                if tries > 100, ph(p) = cand; break; end
            end
        end
    else
        ph(1) = round((-180 + 360*rand)/grid_step)*grid_step;
        for p = 2:P
            while true
                cand = round((-180 + 360*rand)/grid_step)*grid_step;
                if all(abs(wrapTo180(cand - ph(1:p-1))) >= min_sep), ph(p)=cand; break; end
            end
        end
    end
end

function [frob_it, compres_it, ssq_it, c_fin] = multidir_selfcal_one(method, drive_method, ...
    B, Xc_cell, q_cell, phitrue, M, r, lambda, theta_deg, beta_uca, A_dict, phi_grid_deg, ...
    maxIter, C_true, normCtrue, inner_iter, inner_tol, ridge_lambda)
% MULTIDIR_SELFCAL_ONE  Laco alternante multi-dir. As direcoes que ALIMENTAM C
%   vem de 'drive_method' (robusto = KW); a DoA REPORTADA vem de 'method'.
    P = numel(Xc_cell);
    A0 = zeros(M, P);
    for p = 1:P
        phi0 = doa_estimate(Xc_cell{p}, 'KW', q_cell{p}.q, r, lambda, beta_uca, A_dict, phi_grid_deg);
        A0(:,p) = utils.steering_vec_uca(M, r, lambda, theta_deg, phi0);
    end
    C = estimate_C_md(B, A0, M, inner_iter, inner_tol, ridge_lambda);

    frob_it=zeros(maxIter,1); compres_it=zeros(maxIter,1); ssq_it=zeros(maxIter,1);
    same = strcmpi(drive_method, method);
    for it = 1:maxIter
        D = safe_inv(C);
        A = zeros(M, P); ssq = 0;
        for p = 1:P
            Y = D * Xc_cell{p};
            phi_d = doa_estimate(Y, drive_method, q_cell{p}.q, r, lambda, beta_uca, A_dict, phi_grid_deg);
            A(:,p) = utils.steering_vec_uca(M, r, lambda, theta_deg, phi_d);
            if same, phi_e = phi_d;
            else,    phi_e = doa_estimate(Y, method, q_cell{p}.q, r, lambda, beta_uca, A_dict, phi_grid_deg);
            end
            ssq = ssq + wrapTo180(phi_e - phitrue(p))^2;
        end
        C = estimate_C_md(B, A, M, inner_iter, inner_tol, ridge_lambda);
        frob_it(it)    = norm(C - C_true,'fro')/normCtrue;
        compres_it(it) = comp_residual(C, C_true);
        ssq_it(it)     = ssq;
    end
    c_fin = C(1, 2:floor(M/2)+1).';
end

function [C_hat, c_hat, alpha_hat] = estimate_C_md(B_hat, A, M, max_iter, tol, ridge_lambda)
% ESTIMATE_C_MD  Como estimate_C_multidir, mas com Tikhonov (ridge) opcional
%   no LS dos coeficientes (ridge_lambda > 0 doma o mau condicionamento).
    if nargin < 6 || isempty(ridge_lambda), ridge_lambda = 0; end
    [~, P] = size(B_hat);
    K = floor(M/2); is_even = (mod(M,2)==0);
    [C_hat, c_hat, alpha_1, ~] = estimate_C_circulant_uca(B_hat(:,1), A(:,1), M);
    alpha_hat = zeros(P,1); alpha_hat(1) = alpha_1;
    C_prev = C_hat;
    for it = 1:max_iter
        for p = 1:P
            Ca = C_hat * A(:,p);
            alpha_hat(p) = (Ca' * B_hat(:,p)) / (Ca' * Ca);
        end
        M_big = zeros(M*P, K); rhs = zeros(M*P, 1);
        for p = 1:P
            a_p=A(:,p); b_p=B_hat(:,p); ap=alpha_hat(p);
            for i = 1:M
                ridx=(p-1)*M+i;
                for k = 1:K
                    ip=mod(i-1+k,M)+1; im=mod(i-1-k,M)+1;
                    if is_even && (k==K), M_big(ridx,k)=ap*a_p(ip);
                    else,                 M_big(ridx,k)=ap*(a_p(ip)+a_p(im)); end
                end
                rhs(ridx)=b_p(i)-ap*a_p(i);
            end
        end
        if ridge_lambda > 0
            G = M_big'*M_big; f = M_big'*rhs;
            mu = ridge_lambda * trace(G)/K;
            c_hat = (G + mu*eye(K)) \ f;
        else
            c_hat = M_big \ rhs;
        end
        C_new = reconstruct_C(c_hat, M);
        delta = norm(C_new - C_prev,'fro')/max(norm(C_prev,'fro'),eps);
        C_hat = C_new; C_prev = C_new;
        if delta < tol, break; end
    end
end

function C = reconstruct_C(c, M)
% RECONSTRUCT_C  Monta C circulante simetrica a partir dos K coef. unicos.
    K = floor(M/2); is_even = (mod(M,2)==0); c = c(:);
    if is_even, first_row = [1, c(1:K-1).', c(K), flip(c(1:K-1).')];
    else,       first_row = [1, c(1:K).', flip(c(1:K).')]; end
    C = zeros(M);
    for i = 1:M, C(i,:) = circshift(first_row, [0, i-1]); end
end

function rr = comp_residual(C_hat, C_true)
% COMP_RESIDUAL  Residuo de compensacao efetiva: quao perto inv(C_hat)*C_true
%   esta da identidade (apos remover ganho escalar). rr=0 => compensacao perfeita;
%   rr~1 => praticamente sem compensacao.
    M = size(C_hat,1);
    E = safe_inv(C_hat) * C_true;
    g = trace(E)/M;
    if abs(g) > eps, E = E / g; end
    if all(isfinite(E(:)))
        rr = min(norm(E - eye(M),'fro')/sqrt(M), 100);
    else
        rr = 100;
    end
end

function e = scaleinv_frob(C_hat, C_true)
% SCALEINV_FROB  Frobenius relativo apos alinhar o ganho complexo global.
    s = (C_true(:)' * C_hat(:)) / (C_true(:)' * C_true(:));
    e = norm(C_hat - s*C_true,'fro') / norm(C_true,'fro');
end

function D = safe_inv(C)
% SAFE_INV  Inversao robusta a singularidade (nunca retorna Inf/NaN).
    M = size(C,1);
    ws = warning('off','MATLAB:singularMatrix');
    wn = warning('off','MATLAB:nearlySingularMatrix');
    cleanupObj = onCleanup(@() warning([ws wn]));
    rc = rcond(C);
    if isfinite(rc) && rc > 1e-12
        D = inv(C);
        if all(isfinite(D(:))), return; end
    end
    mu = 1e-6 * (norm(C,'fro')^2 / M + eps);
    D  = (C'*C + mu*eye(M)) \ C';
    if ~all(isfinite(D(:))), D = eye(M); end
end

function phi_hat = doa_estimate(X, method, q, r, lambda, beta_uca, A_dict, phi_grid_deg)
% DOA_ESTIMATE  Estima o azimute (graus) por um dos 4 metodos.
    M = size(X,1);
    switch upper(method)
        case 'KW'
            [~, phi] = doa_kw_uca(X, q(:).', r, lambda, beta_uca);
            phi_hat  = phi(1);
        case 'DAS'
            R = (X*X')/size(X,2);
            scores = real(sum(conj(A_dict).*(R*A_dict), 1));
            [~,ip] = max(scores); phi_hat = phi_grid_deg(ip);
        case 'CAPON'
            R = (X*X')/size(X,2);
            R = R + 1e-6*trace(R)/M*eye(M);
            Rinv = R\eye(M);
            den = real(sum(conj(A_dict).*(Rinv*A_dict), 1));
            scores = 1./max(den, eps); [~,ip] = max(scores); phi_hat = phi_grid_deg(ip);
        case 'MUSIC'
            R = (X*X')/size(X,2); R = (R+R')/2;
            [V, Dg] = eig(R); [~, idx] = sort(real(diag(Dg)), 'descend'); V = V(:, idx);
            En = V(:, 2:end);
            proj = En'*A_dict;
            den = real(sum(conj(proj).*proj, 1));
            scores = 1./max(den, eps); [~,ip] = max(scores); phi_hat = phi_grid_deg(ip);
        otherwise
            error('doa_estimate:method', 'Metodo desconhecido: %s', method);
    end
end

function Ctx = compute_Ctx_for_R(fc, M, R, Z0)
% COMPUTE_CTX_FOR_R  Matriz de acoplamento de um UCA de dipolos (parametros-S).
    c = 3e8; lambda = c/fc;
    mp = dipole; mp.Length = 0.5*lambda; mp.Width = 0.01*lambda;
    uca = circularArray; uca.Element = mp; uca.NumElements = M; uca.Radius = R;
    Sobj = sparameters(uca, fc);
    Z_matrix = s2z(Sobj.Parameters(:,:,1), Z0);
    denom = diag(Z_matrix) + Z0;
    Ctx = eye(M);
    for j = 1:M
        for i = 1:M
            if i ~= j, Ctx(i,j) = Z_matrix(i,j) / denom(j); end
        end
    end
end
