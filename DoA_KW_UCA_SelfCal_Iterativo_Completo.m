% =========================================================================
% DoA_KW_UCA_SelfCal_Iterativo_Completo
%
% Single-shot (uma direcao por realizacao) com PASSO DE RELAXACAO 'u'
% ESCOLHIVEL (u=1 => sem suavizacao; u<1 => atenua oscilacoes). Reune TODAS
% as figuras do SelfCal Iterativo e do SelfCal Iterativo Damped, EXCETO a
% varredura de u (efeito_passo_u) -- aqui u e' fixo, definido no topo.
%
%   Laco alternante (init = KW):  C_t = (1-u)*C_{t-1} + u*C_LS
%   4 metodos de DoA em paralelo (KW/DAS/Capon/MUSIC).
%   ORACLE iterativo (ang. verdadeiro) como referencia/limite.
%
%   FIGURAS:
%     (A) convergencia por SNR: DoA RMSE, Frobenius e residuo de compensacao
%         vs iteracao (com oracle).
%     (B) SCATTER de erro por azimute (iterativo vs sem comp.), 1 painel/metodo.
%     (C) resumo de DoA vs SNR: RMSE e MEDIANA (iterativo vs sem comp.).
%     (D) resumo de metricas de C vs SNR: Frobenius bruto, invariante a escala
%         e residuo de compensacao (com oracle).
%     (E) PLANO COMPLEXO dos coeficientes (verdadeiro vs nuvem estimada + oracle).
% =========================================================================

clear; clc; close all;

% Avisos ESPERADOS em abertura pequena (steering quase paralelos -> LS de C
% mal-condicionado, alpha~0 -> coeficientes nao-finitos). Sao tratados: as
% estimativas nao-finitas sao DESCARTADAS das medias. Suprimimos o ruido.
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:illConditionedMatrix');
warning('off','estimate_C_circulant_uca:smallAlpha');

%% ---- Parametros do array ----
M      = 8;  fc = 500e6;  c = 3e8;  lambda = c/fc;
r      = 0.15*lambda;      % <-- raio (escolhivel)
theta_sig_deg = 90;

%% ---- Parametros do sinal (FSK-2 conhecido) ----
fs = 288000; N = 9900; Rs = 9600; sps = 30; alpha = 0.3; span = 8; fd = 4.8e3;

%% ---- Parametros do experimento ----
maxIter      = 6;
u            = 1;                 % <-- PASSO de relaxacao ESCOLHIVEL (1 = sem suavizar)
range_SNR_dB = [-6, 0, 6, 12];
methods      = {'KW','DAS','CAPON','MUSIC'};
nMethods     = numel(methods);
phi_grid_deg = -180:0.1:180;
grid_step    = phi_grid_deg(2) - phi_grid_deg(1);
beta_uca     = 2*pi*(0:M-1).'/M;
K            = floor(M/2);

% --- Monte Carlo ---
n_angles = 100;             % azimutes aleatorios (no grid)
n_trials = 5;             % realizacoes por azimute
n_real   = n_angles*n_trials;
rng(2026, 'twister');
phi_set  = round( (-180 + 360*rand(1, n_angles)) / grid_step ) * grid_step;

outDir = fullfile(pwd, 'Iterative Completo Graphs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---- Matriz de acoplamento VERDADEIRA (so' p/ medir erro) ----
C_true     = compute_Ctx_for_R(fc, M, r, 50);
normCtrue  = norm(C_true, 'fro');
c_true_vec = C_true(1, 2:K+1).';
fprintf('C_true (raio=%.2f lambda, ||C_true||_F=%.4f, cond=%.1f),  u=%.2f\n', ...
    r/lambda, normCtrue, cond(C_true), u);

%% ---- Dicionario de steering vectors ----
nGrid = numel(phi_grid_deg);
A_dict = zeros(M, nGrid);
for ig = 1:nGrid
    A_dict(:,ig) = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_grid_deg(ig));
end

%% ---- Acumuladores ----
nSNR             = numel(range_SNR_dB);
sumsq_doa_iter   = zeros(nMethods, maxIter, nSNR);   % DoA RMSE por iteracao
sum_frob_iter    = zeros(nMethods, maxIter, nSNR);   % Frobenius por iteracao
sum_cres_iter    = zeros(nMethods, maxIter, nSNR);   % residuo de compensacao por iteracao
sum_frob_oracle  = zeros(maxIter, nSNR);             % oracle (trajetoria)
sum_cres_oracle  = zeros(maxIter, nSNR);
c_runs           = zeros(nMethods, K, n_real, nSNR); % coef. finais (plano complexo)
c_oracle         = zeros(K, n_real, nSNR);
Eabs_final  = zeros(nMethods, n_angles, n_trials, nSNR);  % |erro DoA| ultima iter (scatter/mediana)
Eabs_noComp = zeros(nMethods, n_angles, n_trials, nSNR);  % |erro DoA| sem compensacao
Eabs_noMC   = zeros(nMethods, n_angles, n_trials, nSNR);  % |erro DoA| sem acoplamento
cnt         = zeros(1, nSNR);
cnt_fi      = zeros(nMethods, maxIter, nSNR);   % nº de validos (frob/cres por iteracao, metodos)
cnt_fo      = zeros(maxIter, nSNR);             % nº de validos do oracle

%% ======================= LOOP PRINCIPAL =======================
t0 = tic;
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fprintf('\n===== SNR = %+d dB  (u=%.2f) =====\n', SNR_dB, u);
    ir = 0;
    for ia = 1:n_angles
        phi_true = phi_set(ia);
        for itr = 1:n_trials
            ir = ir + 1;
            [~, q, ~, Xsig, ~, Xn] = utils.simulate_fsk_data_uca(M, r, lambda, ...
                phi_true, phi_true+90, theta_sig_deg, theta_sig_deg, ...
                SNR_dB, -100, N, fs, Rs, sps, alpha, span, fd);
            q = q(:);
            X_ideal    = Xsig + Xn;
            X_coupled  = C_true*Xsig + Xn;
            b_hat_orig = X_coupled*conj(q)/(q'*q);

            % --- DoA de referencia (sem acoplamento e sem compensacao) ---
            for im = 1:nMethods
                phiMC = doa_estimate(X_ideal,   methods{im}, q, r, lambda, beta_uca, A_dict, phi_grid_deg);
                phiNC = doa_estimate(X_coupled, methods{im}, q, r, lambda, beta_uca, A_dict, phi_grid_deg);
                Eabs_noMC(im,ia,itr,iSNR)   = abs(wrapTo180(phiMC - phi_true));
                Eabs_noComp(im,ia,itr,iSNR) = abs(wrapTo180(phiNC - phi_true));
            end

            % --- ORACLE iterativo (angulo verdadeiro, mesmo laco/u) ---
            [frob_orc_it, cres_orc_it, ~, c_orc_fin] = selfcal_one_damped('ORACLE', ...
                X_coupled, b_hat_orig, q, phi_true, M, r, lambda, theta_sig_deg, ...
                beta_uca, A_dict, phi_grid_deg, C_true, normCtrue, u, maxIter);
            fo = frob_orc_it(:); co = cres_orc_it(:);
            vfo = isfinite(fo) & isfinite(co);  fo(~vfo)=0; co(~vfo)=0;   % descarta nao-finitos
            sum_frob_oracle(:,iSNR) = sum_frob_oracle(:,iSNR) + fo;
            sum_cres_oracle(:,iSNR) = sum_cres_oracle(:,iSNR) + co;
            cnt_fo(:,iSNR)          = cnt_fo(:,iSNR) + vfo;
            c_oracle(:,ir,iSNR)     = c_orc_fin(:);

            % --- Self-cal DAMPED, por metodo ---
            for im = 1:nMethods
                [frob_it, cres_it, derr_it, c_fin] = selfcal_one_damped(methods{im}, ...
                    X_coupled, b_hat_orig, q, phi_true, M, r, lambda, theta_sig_deg, ...
                    beta_uca, A_dict, phi_grid_deg, C_true, normCtrue, u, maxIter);
                fi = frob_it(:).'; ci = cres_it(:).';
                vfi = isfinite(fi) & isfinite(ci);  fi(~vfi)=0; ci(~vfi)=0;   % descarta nao-finitos
                sum_frob_iter(im,:,iSNR)  = sum_frob_iter(im,:,iSNR)  + fi;
                sum_cres_iter(im,:,iSNR)  = sum_cres_iter(im,:,iSNR)  + ci;
                cnt_fi(im,:,iSNR)         = cnt_fi(im,:,iSNR) + vfi;
                sumsq_doa_iter(im,:,iSNR) = sumsq_doa_iter(im,:,iSNR) + (derr_it(:).').^2;
                Eabs_final(im,ia,itr,iSNR)= derr_it(end);
                c_runs(im,:,ir,iSNR)      = c_fin(:).';
            end
            cnt(iSNR) = cnt(iSNR) + 1;
        end
        fprintf('  phi=%+6.1f deg  | decorrido %.0fs\n', phi_true, toc(t0));
    end
end

%% ---- Consolidacao ----
rmse_doa_iter  = zeros(nMethods, maxIter, nSNR);
mean_frob_iter = zeros(nMethods, maxIter, nSNR);
mean_cres_iter = zeros(nMethods, maxIter, nSNR);
for iSNR = 1:nSNR
    rmse_doa_iter(:,:,iSNR)  = sqrt(sumsq_doa_iter(:,:,iSNR) / cnt(iSNR));   % DoA sempre finito
    mean_frob_iter(:,:,iSNR) = sum_frob_iter(:,:,iSNR) ./ max(cnt_fi(:,:,iSNR),1);   % media robusta
    mean_cres_iter(:,:,iSNR) = sum_cres_iter(:,:,iSNR) ./ max(cnt_fi(:,:,iSNR),1);
end
rmse_doa_final   = squeeze(rmse_doa_iter(:,maxIter,:));
frob_final       = squeeze(mean_frob_iter(:,maxIter,:));
cres_final       = squeeze(mean_cres_iter(:,maxIter,:));
mean_frob_oracle = sum_frob_oracle ./ max(cnt_fo,1);   % (maxIter x nSNR) robusto
mean_cres_oracle = sum_cres_oracle ./ max(cnt_fo,1);
frob_oracle      = mean_frob_oracle(end,:);
cres_oracle      = mean_cres_oracle(end,:);

% Por azimute (RMSE/mediana sobre trials) e global
rmse_angle_iter   = squeeze(sqrt(mean(Eabs_final.^2,  3)));   % nMethods x n_angles x nSNR
rmse_angle_noComp = squeeze(sqrt(mean(Eabs_noComp.^2, 3)));
E2_final  = reshape(Eabs_final,  nMethods, n_real, nSNR);
E2_noComp = reshape(Eabs_noComp, nMethods, n_real, nSNR);
E2_noMC   = reshape(Eabs_noMC,   nMethods, n_real, nSNR);
rmse_doa_noComp = squeeze(sqrt(mean(E2_noComp.^2, 2)));
rmse_doa_noMC   = squeeze(sqrt(mean(E2_noMC.^2,   2)));
med_final       = squeeze(median(E2_final,  2));
med_noComp      = squeeze(median(E2_noComp, 2));

% Metricas de C derivadas dos coeficientes
scaleinv_final = zeros(nMethods, nSNR);  coeff_err     = zeros(nMethods, K, nSNR);
scaleinv_orc   = zeros(1, nSNR);         coeff_err_orc = zeros(K, nSNR);
% Medias ROBUSTAS (ignoram realizacoes com coeficientes nao-finitos -> sem NaN)
for iSNR = 1:nSNR
    valsO = zeros(n_real,1);
    for run = 1:n_real
        valsO(run) = scaleinv_frob(reconstruct_C(c_oracle(:,run,iSNR), M), C_true);
    end
    scaleinv_orc(iSNR) = rmean(valsO);
    for k = 1:K, coeff_err_orc(k,iSNR) = rmean(abs(c_oracle(k,:,iSNR).' - c_true_vec(k))); end
    for im = 1:nMethods
        vals = zeros(n_real,1);
        for run = 1:n_real
            vals(run) = scaleinv_frob(reconstruct_C(squeeze(c_runs(im,:,run,iSNR)).', M), C_true);
        end
        scaleinv_final(im,iSNR) = rmean(vals);
        for k = 1:K, coeff_err(im,k,iSNR) = rmean(abs(squeeze(c_runs(im,k,:,iSNR)) - c_true_vec(k))); end
    end
end

flr = @(x) max(x, 1e-3);
colorsM   = lines(nMethods);
markers_m = {'o-','s-','d-','^-'};
mk_only   = {'o','s','d','^'};
ck_colors = [0.85 0.10 0.10; 0.10 0.60 0.15; 0.10 0.30 0.85; 0.60 0.15 0.65];

%% ---- (A) convergencia por SNR: DoA, Frobenius, residuo de compensacao ----
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fig = figure('Name',sprintf('Convergencia SNR=%+d dB',SNR_dB),'Color','w','Position',[40 60 1500 470]);

    subplot(1,3,1); hold on; grid on;
    for im = 1:nMethods
        plot(1:maxIter, flr(squeeze(rmse_doa_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:),'LineWidth',1.6,'MarkerFaceColor',colorsM(im,:),'DisplayName',methods{im});
        yline(flr(rmse_doa_noComp(im,iSNR)), ':','Color',colorsM(im,:),'LineWidth',1.0,'HandleVisibility','off');
    end
    set(gca,'YScale','log'); xlabel('Iteracao'); ylabel('RMSE de DoA (graus)');
    title('DoA (pontilhado = sem comp.)'); legend('Location','best'); xlim([1 maxIter]);

    subplot(1,3,2); hold on; grid on;
    for im = 1:nMethods
        plot(1:maxIter, flr(squeeze(mean_frob_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:),'LineWidth',1.6,'MarkerFaceColor',colorsM(im,:),'DisplayName',methods{im});
    end
    plot(1:maxIter, flr(mean_frob_oracle(:,iSNR)), 'k--p','LineWidth',1.6,'MarkerFaceColor','k','MarkerSize',5,'DisplayName','oracle (ang. sabido)');
    set(gca,'YScale','log'); xlabel('Iteracao'); ylabel('Frobenius ||C_{hat}-C_{true}||/||C_{true}||');
    title('Frobenius de C vs iteracao'); legend('Location','best'); xlim([1 maxIter]);

    subplot(1,3,3); hold on; grid on;
    for im = 1:nMethods
        plot(1:maxIter, flr(squeeze(mean_cres_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:),'LineWidth',1.6,'MarkerFaceColor',colorsM(im,:),'DisplayName',methods{im});
    end
    plot(1:maxIter, flr(mean_cres_oracle(:,iSNR)), 'k--p','LineWidth',1.6,'MarkerFaceColor','k','MarkerSize',5,'DisplayName','oracle (ang. sabido)');
    set(gca,'YScale','log'); xlabel('Iteracao'); ylabel('residuo de compensacao');
    title('Residuo de compensacao vs iteracao'); legend('Location','best'); xlim([1 maxIter]);

    sgtitle(sprintf('Convergencia  SNR=%+d dB, u=%.2f, raio=%.2f\\lambda', SNR_dB, u, r/lambda),'FontWeight','bold');
    exportgraphics(fig, fullfile(outDir, sprintf('convergencia_SNR_%+03d.png',SNR_dB)),'Resolution',170);
end

%% ---- (B) SCATTER de erro por azimute (1 painel/metodo, mesmo ylim) ----
[phi_sorted, ord] = sort(phi_set);
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fig = figure('Name',sprintf('Erro por azimute SNR=%+d dB',SNR_dB),'Color','w','Position',[60 60 1300 820]);
    vv   = flr([reshape(rmse_angle_iter(:,:,iSNR),[],1); reshape(rmse_angle_noComp(:,:,iSNR),[],1)]);
    yl_s = [min(vv)*0.8, max(vv)*1.3];
    for im = 1:nMethods
        subplot(2,2,im); hold on; grid on;
        semilogy(phi_sorted, flr(rmse_angle_iter(im,ord,iSNR)), [mk_only{im} '-'], ...
            'Color',colorsM(im,:),'LineWidth',1.4,'MarkerFaceColor',colorsM(im,:),'DisplayName','iterativo');
        semilogy(phi_sorted, flr(rmse_angle_noComp(im,ord,iSNR)), [mk_only{im} ':'], ...
            'Color',colorsM(im,:),'LineWidth',1.0,'MarkerFaceColor','none','DisplayName','sem comp.');
        yline(flr(med_final(im,iSNR)), 'k--', 'LineWidth',1.0, 'DisplayName','mediana global (iter)');
        set(gca,'YScale','log'); xlabel('azimute \phi (graus)'); ylabel('RMSE de DoA (graus)');
        xlim([-180 180]); xticks(-180:90:180); ylim(yl_s);
        title(sprintf('%s  | mediana iter=%.3f, semComp=%.3f', methods{im}, med_final(im,iSNR), med_noComp(im,iSNR)));
        legend('Location','best');
    end
    sgtitle(sprintf('Erro de DoA por azimute  (SNR=%+d dB, u=%.2f, %d realizacoes/azimute)', SNR_dB, u, n_trials),'FontWeight','bold');
    exportgraphics(fig, fullfile(outDir, sprintf('scatter_azimute_SNR_%+03d.png',SNR_dB)),'Resolution',170);
end

%% ---- (C) resumo de DoA vs SNR: RMSE e MEDIANA ----
plot_summary_vs_snr(range_SNR_dB, rmse_doa_final, rmse_doa_noComp, methods, colorsM, markers_m, flr, ...
    n_angles, n_trials, r/lambda, u, 'RMSE de DoA', fullfile(outDir,'resumo_RMSE_doa_vs_snr.png'));
plot_summary_vs_snr(range_SNR_dB, med_final, med_noComp, methods, colorsM, markers_m, flr, ...
    n_angles, n_trials, r/lambda, u, 'MEDIANA de erro de DoA', fullfile(outDir,'resumo_MEDIANA_doa_vs_snr.png'));

%% ---- (D) resumo de metricas de C vs SNR ----
fig = figure('Color','w','Position',[60 80 1500 470]);
metsD = {frob_final, scaleinv_final, cres_final};
orcD  = {frob_oracle, scaleinv_orc, cres_oracle};
titD  = {'Frobenius bruto', 'Frobenius invariante a escala', 'Residuo de compensacao'};
for sp = 1:3
    subplot(1,3,sp); hold on; grid on;
    for im = 1:nMethods
        plot(range_SNR_dB, flr(metsD{sp}(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
            'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',methods{im});
    end
    plot(range_SNR_dB, flr(orcD{sp}), 'k--p', 'LineWidth',1.6,'MarkerFaceColor','k','MarkerSize',8,'DisplayName','oracle (ang. sabido)');
    set(gca,'YScale','log'); xlabel('SNR (dB)'); ylabel('erro de C'); xticks(range_SNR_dB);
    title(titD{sp}); legend('Location','best');
end
sgtitle(sprintf('Metricas de C vs SNR  (u=%.2f, %d realizacoes, raio=%.2f\\lambda)', u, n_real, r/lambda),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'resumo_metricas_C_vs_snr.png'),'Resolution',170);

%% ---- (E) plano complexo dos coeficientes por SNR ----
for iSNR = 1:nSNR
    fig = figure('Name',sprintf('Coef. complexos SNR=%+d',range_SNR_dB(iSNR)),'Color','w','Position',[40 40 1500 820]);
    srcC = cell(1, nMethods+1);
    for im = 1:nMethods, srcC{im} = squeeze(c_runs(im,:,:,iSNR)); end
    srcC{nMethods+1} = squeeze(c_oracle(:,:,iSNR));
    labC = [methods, {'ORACLE (ang. sabido)'}];
    errC = [coeff_err(:,:,iSNR); coeff_err_orc(:,iSNR).'];
    lim  = common_complex_lim(srcC, c_true_vec);
    for pp = 1:(nMethods+1)
        subplot(2,3,pp); hold on; grid on; axis equal;
        S = srcC{pp};
        for k = 1:K
            scatter(real(S(k,:)), imag(S(k,:)), 16, ck_colors(k,:), 'filled', ...
                'MarkerFaceAlpha',0.30, 'DisplayName',sprintf('c_%d (est)',k));
        end
        for k = 1:K
            plot(real(c_true_vec(k)), imag(c_true_vec(k)), 'p', 'MarkerSize',16, ...
                'MarkerFaceColor',ck_colors(k,:),'MarkerEdgeColor','k','LineWidth',1.3,'HandleVisibility','off');
            text(real(c_true_vec(k)), imag(c_true_vec(k)), sprintf(' c_%d',k),'FontWeight','bold','FontSize',10);
        end
        xlabel('Re'); ylabel('Im'); xlim(lim); ylim(lim);
        estr = strjoin(arrayfun(@(k) sprintf('%.2f',errC(pp,k)), 1:K, 'uni',0), ', ');
        title(sprintf('%s  |\\Deltac_k|=[%s]', labC{pp}, estr), 'FontSize',9);
        if pp == 1, legend('Location','bestoutside'); end
    end
    sgtitle(sprintf(['Coeficientes c_1..c_%d: verdadeiro (estrela) vs estimado (nuvem)\n' ...
        'SNR=%+d dB, u=%.2f, %d realizacoes, raio=%.2f\\lambda'], ...
        K, range_SNR_dB(iSNR), u, n_real, r/lambda),'FontWeight','bold');
    exportgraphics(fig, fullfile(outDir, sprintf('coeficientes_complexo_SNR_%+03d.png',range_SNR_dB(iSNR))),'Resolution',160);
end

%% ---- (F) BARRAS: |Delta c_k| por coeficiente, comparando metodos (+oracle) ----
% Gráfico com VALORES (alturas das barras) para comparar os metodos por
% coeficiente, complementando a nuvem do plano complexo. 1 painel por SNR.
labB = [methods, {'oracle'}];
allb = [coeff_err(:); coeff_err_orc(:)];  allb = allb(isfinite(allb) & allb>0);
if isempty(allb), ylb = [1e-3 1]; else, ylb = [min(allb)*0.5, max(allb)*1.5]; end
fig = figure('Name','Erro por coeficiente (barras)','Color','w','Position',[60 60 1300 820]);
for iSNR = 1:nSNR
    subplot(2,2,iSNR);
    barData = [coeff_err(:,:,iSNR); coeff_err_orc(:,iSNR).'];   % (nMethods+1) x K
    hb = bar(1:K, barData.');                                  % grupos=coef, barras=metodo+oracle
    set(hb, 'BaseValue', ylb(1));  grid on; set(gca,'YScale','log'); ylim(ylb);
    xticks(1:K); xticklabels(arrayfun(@(k)sprintf('c_%d',k),1:K,'uni',0));
    xlabel('coeficiente'); ylabel('|\Deltac_k| medio');
    title(sprintf('SNR = %+d dB', range_SNR_dB(iSNR)));
    if iSNR==1, legend(hb, labB, 'Location','best'); end
end
sgtitle(sprintf('Erro medio por coeficiente |\\Deltac_k| (valores)  (u=%.2f, raio=%.2f\\lambda)', ...
    u, r/lambda),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'erro_por_coeficiente_barras.png'),'Resolution',170);

%% ---- Resumo numerico ----
fprintf('\n===== RESUMO (u=%.2f, %d realizacoes) =====\n', u, n_real);
for iSNR = 1:nSNR
    fprintf('SNR=%+3d dB:   [oracle] Frob=%.3f FrobEsc=%.3f compRes=%.3f |dc_k|=[%s]\n', ...
        range_SNR_dB(iSNR), frob_oracle(iSNR), scaleinv_orc(iSNR), cres_oracle(iSNR), ...
        strjoin(arrayfun(@(k) sprintf('%.2f',coeff_err_orc(k,iSNR)),1:K,'uni',0),', '));
    for im = 1:nMethods
        fprintf('   %-6s  Frob=%.3f FrobEsc=%.3f compRes=%.3f  |dc_k|=[%s]  DoA(rmse/med)=%.3f/%.3f\n', ...
            methods{im}, frob_final(im,iSNR), scaleinv_final(im,iSNR), cres_final(im,iSNR), ...
            strjoin(arrayfun(@(k) sprintf('%.2f',coeff_err(im,k,iSNR)),1:K,'uni',0),', '), ...
            rmse_doa_final(im,iSNR), med_final(im,iSNR));
    end
end
fprintf('\nFiguras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function plot_summary_vs_snr(range_SNR_dB, M_iter, M_noComp, methods, colorsM, markers_m, flr, ...
    n_angles, n_trials, r_over_lambda, u, ylab, fname)
    nMethods = numel(methods);
    fig = figure('Color','w','Position',[100 100 920 560]); hold on; grid on;
    for im = 1:nMethods
        plot(range_SNR_dB, flr(M_iter(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
            'LineWidth',2.0,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',sprintf('%s iterativo',methods{im}));
        plot(range_SNR_dB, flr(M_noComp(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
            'LineWidth',1.2,'LineStyle',':','MarkerFaceColor','none','MarkerSize',8,'DisplayName',sprintf('%s sem comp.',methods{im}));
    end
    set(gca,'YScale','log'); xlabel('SNR (dB)'); ylabel(sprintf('%s (graus)', ylab));
    xticks(range_SNR_dB); xlim([min(range_SNR_dB)-1, max(range_SNR_dB)+1]);
    title(sprintf(['%s: laco iterativo (solido) vs sem compensacao (pontilhado)\n' ...
                   '%d azimutes x %d realizacoes,  u=%.2f,  raio=%.2f\\lambda'], ...
                   ylab, n_angles, n_trials, u, r_over_lambda));
    legend('Location','eastoutside');
    exportgraphics(fig, fname, 'Resolution',180);
end

function [frob_it, cres_it, derr_it, c_fin] = selfcal_one_damped(method, X_coupled, ...
    b_hat_orig, q, phi_true, M, r, lambda, theta_deg, beta_uca, A_dict, phi_grid_deg, ...
    C_true, normCtrue, u, maxIter)
% Laco alternante single-shot com passo de relaxacao u. method='ORACLE' usa
% sempre o angulo verdadeiro (DoA perfeita) no lugar da DoA estimada.
    phi0  = doa_estimate(X_coupled, 'KW', q, r, lambda, beta_uca, A_dict, phi_grid_deg);
    a0    = utils.steering_vec_uca(M, r, lambda, theta_deg, phi0);
    C_hat = estimate_C_circulant_uca(b_hat_orig, a0, M);
    frob_it=zeros(maxIter,1); cres_it=zeros(maxIter,1); derr_it=zeros(maxIter,1);
    for it = 1:maxIter
        D   = safe_inv(C_hat);
        Y   = D * X_coupled;
        if strcmpi(method,'ORACLE')
            phi = phi_true;
        else
            phi = doa_estimate(Y, method, q, r, lambda, beta_uca, A_dict, phi_grid_deg);
        end
        a_hat = utils.steering_vec_uca(M, r, lambda, theta_deg, phi);
        C_ls  = estimate_C_circulant_uca(b_hat_orig, a_hat, M);
        C_hat = (1-u)*C_hat + u*C_ls;
        frob_it(it) = norm(C_hat - C_true,'fro')/normCtrue;
        cres_it(it) = comp_residual(C_hat, C_true);
        derr_it(it) = abs(wrapTo180(phi - phi_true));
    end
    c_fin = C_hat(1, 2:floor(M/2)+1).';
end

function lim = common_complex_lim(srcC, c_true_vec)
    parts = cell(numel(srcC)+1,1);
    for pp = 1:numel(srcC)
        S = srcC{pp};  parts{pp} = [real(S(:)); imag(S(:))];
    end
    parts{end} = [real(c_true_vec(:)); imag(c_true_vec(:))];
    v = vertcat(parts{:});  v = sort(v(isfinite(v)));  n = numel(v);
    if n < 5, lim = [-1 1]; return; end
    lo = v(max(1,round(0.02*n)));  hi = v(min(n,round(0.98*n)));
    lo = min(lo, min(parts{end}));  hi = max(hi, max(parts{end}));
    pad = 0.1*(hi-lo) + eps;  lim = [lo-pad, hi+pad];
end

function m = rmean(x)
% Media ROBUSTA: ignora elementos nao-finitos (NaN/Inf). NaN so' se TUDO falhar.
    x = x(isfinite(x));
    if isempty(x), m = NaN; else, m = mean(x); end
end

function C = reconstruct_C(c, M)
    K = floor(M/2); is_even = (mod(M,2)==0); c = c(:);
    if is_even, first_row = [1, c(1:K-1).', c(K), flip(c(1:K-1).')];
    else,       first_row = [1, c(1:K).', flip(c(1:K).')]; end
    C = zeros(M);
    for i = 1:M, C(i,:) = circshift(first_row, [0, i-1]); end
end

function rr = comp_residual(C_hat, C_true)
    if ~all(isfinite(C_hat(:))), rr = NaN; return; end
    M = size(C_hat,1);
    E = safe_inv(C_hat) * C_true;
    g = trace(E)/M; if abs(g) > eps, E = E / g; end
    if all(isfinite(E(:))), rr = norm(E - eye(M),'fro')/sqrt(M);
    else,                   rr = NaN; end
end

function e = scaleinv_frob(C_hat, C_true)
    s = (C_true(:)' * C_hat(:)) / (C_true(:)' * C_true(:));
    e = norm(C_hat - s*C_true,'fro') / norm(C_true,'fro');
end

function D = safe_inv(C)
    M = size(C,1);
    ws = warning('off','MATLAB:singularMatrix');
    wn = warning('off','MATLAB:nearlySingularMatrix');
    cleanupObj = onCleanup(@() warning([ws wn]));
    rc = rcond(C);
    if isfinite(rc) && rc > 1e-12
        D = inv(C); if all(isfinite(D(:))), return; end
    end
    mu = 1e-6 * (norm(C,'fro')^2 / M + eps);
    D  = (C'*C + mu*eye(M)) \ C';
    if ~all(isfinite(D(:))), D = eye(M); end
end

function phi_hat = doa_estimate(X, method, q, r, lambda, beta_uca, A_dict, phi_grid_deg)
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
            R = (X*X')/size(X,2); R = R + 1e-6*trace(R)/M*eye(M);
            Rinv = R\eye(M);
            den = real(sum(conj(A_dict).*(Rinv*A_dict), 1));
            scores = 1./max(den, eps); [~,ip] = max(scores); phi_hat = phi_grid_deg(ip);
        case 'MUSIC'
            R = (X*X')/size(X,2); R = (R+R')/2;
            [V, Dg] = eig(R); [~, idx] = sort(real(diag(Dg)), 'descend'); V = V(:, idx);
            En = V(:, 2:end); proj = En'*A_dict;
            den = real(sum(conj(proj).*proj, 1));
            scores = 1./max(den, eps); [~,ip] = max(scores); phi_hat = phi_grid_deg(ip);
        otherwise
            error('doa_estimate:method', 'Metodo desconhecido: %s', method);
    end
end

function Ctx = compute_Ctx_for_R(fc, M, R, Z0)
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
