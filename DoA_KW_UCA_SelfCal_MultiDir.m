% =========================================================================
% DoA_KW_UCA_SelfCal_MultiDir
%
% Self-calibration iterativa por COMPENSACAO de acoplamento usando MULTIPLAS
% DIRECOES NAO-SIMULTANEAS. Motivacao: no caso single-shot (uma unica
% direcao por estimativa) o problema "estimar C circulante a partir de uma
% fonte conhecida" e' fracamente identificavel, e o laco alternante so' fica
% estavel quando dirigido pelo KW. Com P direcoes DISTINTAS (medidas em
% momentos diferentes, fonte uma de cada vez), a estimacao conjunta de C fica
% sobredeterminada e deve estabilizar tambem DAS/Capon/MUSIC.
%
%   Modelo (P medidas, nao-simultaneas):
%       b_p = alpha_p * C * a(phi_p) + ruido,   p = 1..P    (C circulante)
%
%   Laco alternante (init = KW, sem damping/parada/oraculo):
%     1) para cada p: phi_p = DoA( inv(C) * X_p )   (4 metodos em paralelo)
%     2) C <- estimate_C_multidir( {b_p}, {a(phi_p)} )  (LS conjunto)
%     3) repete
%
%   MONTE CARLO (sugestao 3): MUITAS rodadas aleatorias, cada uma com P
%   azimutes sorteados (no grid, com separacao minima) e ruido novo.
%
%   Figuras:
%     (A) convergencia por SNR: Frobenius e RMSE de DoA vs iteracao
%     (B) varredura do nº de direcoes P: erro final vs P (mostra o ganho
%         de identificabilidade de usar mais direcoes)
%     (C) resumo vs SNR: multi-dir (solido) vs sem compensacao (pontilhado)
% =========================================================================

clear; clc; close all;

%% ---- Parametros do array ----
M      = 8;            % nº de elementos do UCA
fc     = 500e6;        % Hz
c      = 3e8;
lambda = c/fc;
r      = 0.1*lambda;   % raio 0.2 lambda
theta_sig_deg = 90;    % elevacao (plano XY)

%% ---- Parametros do sinal (FSK-2 conhecido) ----
fs = 288000; N = 2100; Rs = 9600; sps = 30; alpha = 0.3; span = 8; fd = 4.8e3;

%% ---- Parametros do experimento ----
maxIter      = 12;                  % <-- LIMITE de iteracoes do laco externo
range_SNR_dB = [-6, 0, 6, 12];      % cenarios de SNR
methods      = {'KW','DAS','CAPON','MUSIC'};
nMethods     = numel(methods);
phi_grid_deg = -180:0.5:180;
grid_step    = phi_grid_deg(2) - phi_grid_deg(1);
beta_uca     = 2*pi*(0:M-1).'/M;

% --- Multi-direcao ---
P_default    = 5;            % nº de direcoes (medidas nao-simultaneas) por rodada
range_P      = [1 2 3 4 6 8];% varredura do nº de direcoes (figura B)
snr_Psweep   = 6;           % SNR usada na varredura de P (dB)
min_sep_deg  = 10;          % separacao angular minima entre as P direcoes
inner_iter   = 30;          % iteracoes internas do estimate_C_multidir
inner_tol    = 1e-10;

% --- Monte Carlo ---  (sugestao 3: MUITAS rodadas aleatorias, parametrizado)
n_runs   = 60;              % <-- nº de rodadas (cada uma: P direcoes novas + ruido)
rng(2026, 'twister');       % reprodutibilidade

outDir = fullfile(pwd, 'MultiDir SelfCal Graphs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---- Matriz de acoplamento VERDADEIRA (so' p/ medir erro) ----
C_true    = compute_Ctx_for_R(fc, M, r, 50);
normCtrue = norm(C_true, 'fro');
fprintf('C_true (raio=%.2f lambda, ||C_true||_F=%.4f)\n', r/lambda, normCtrue);

%% ---- Dicionario de steering vectors ----
nGrid = numel(phi_grid_deg);
A_dict = zeros(M, nGrid);
for ig = 1:nGrid
    A_dict(:,ig) = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_grid_deg(ig));
end

%% =====================================================================
%   PARTE A/C: convergencia por SNR e resumo vs SNR  (P = P_default)
%% =====================================================================
nSNR            = numel(range_SNR_dB);
sumsq_doa_iter  = zeros(nMethods, maxIter, nSNR);  % SOMA (erro DoA)^2 por iteracao (sobre rodadas x shots)
sum_frob_iter   = zeros(nMethods, maxIter, nSNR);  % SOMA Frobenius por iteracao (sobre rodadas)
sumsq_doa_noComp= zeros(nMethods, nSNR);           % referencia sem compensacao
sumsq_doa_noMC  = zeros(nMethods, nSNR);           % referencia sem acoplamento
cnt_doa         = zeros(1, nSNR);                  % nº de (rodada x shot)
cnt_run         = zeros(1, nSNR);                  % nº de rodadas

t0 = tic;
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fprintf('\n===== [A/C] SNR = %+d dB  (P=%d, %d rodadas) =====\n', SNR_dB, P_default, n_runs);
    for run = 1:n_runs
        % --- sorteia P direcoes distintas e gera os P shots (nao-simultaneos) ---
        [B, Xc_cell, q_cell, phitrue] = gen_shots(P_default, grid_step, min_sep_deg, ...
            M, r, lambda, theta_sig_deg, SNR_dB, N, fs, Rs, sps, alpha, span, fd, C_true);

        % --- referencias (sem iterar): sem acoplamento e sem compensacao ---
        for im = 1:nMethods
            for p = 1:P_default
                X_ideal = q_cell{p}.X_ideal;
                phiMC = doa_estimate(X_ideal,    methods{im}, q_cell{p}.q, r, lambda, beta_uca, A_dict, phi_grid_deg);
                phiNC = doa_estimate(Xc_cell{p}, methods{im}, q_cell{p}.q, r, lambda, beta_uca, A_dict, phi_grid_deg);
                sumsq_doa_noMC(im,iSNR)   = sumsq_doa_noMC(im,iSNR)   + wrapTo180(phiMC - phitrue(p))^2;
                sumsq_doa_noComp(im,iSNR) = sumsq_doa_noComp(im,iSNR) + wrapTo180(phiNC - phitrue(p))^2;
            end
        end

        % --- laco multi-dir por metodo ---
        for im = 1:nMethods
            [frob_it, ssq_it] = multidir_selfcal_one(methods{im}, B, Xc_cell, q_cell, ...
                phitrue, M, r, lambda, theta_sig_deg, beta_uca, A_dict, phi_grid_deg, ...
                maxIter, C_true, normCtrue, inner_iter, inner_tol);
            sum_frob_iter(im,:,iSNR)  = sum_frob_iter(im,:,iSNR)  + frob_it(:).';
            sumsq_doa_iter(im,:,iSNR) = sumsq_doa_iter(im,:,iSNR) + ssq_it(:).';
        end
        cnt_doa(iSNR) = cnt_doa(iSNR) + P_default;
        cnt_run(iSNR) = cnt_run(iSNR) + 1;
    end
    fprintf('  decorrido %.0fs\n', toc(t0));
end

% --- Consolidacao A/C ---
rmse_doa_iter = zeros(nMethods, maxIter, nSNR);
mean_frob_iter= zeros(nMethods, maxIter, nSNR);
for iSNR = 1:nSNR
    rmse_doa_iter(:,:,iSNR)  = sqrt(sumsq_doa_iter(:,:,iSNR) / cnt_doa(iSNR));
    mean_frob_iter(:,:,iSNR) =      sum_frob_iter(:,:,iSNR)  / cnt_run(iSNR);
end
rmse_doa_final = squeeze(rmse_doa_iter(:,maxIter,:));            % nMethods x nSNR
rmse_doa_noComp= sqrt(sumsq_doa_noComp ./ cnt_doa);             % nMethods x nSNR
rmse_doa_noMC  = sqrt(sumsq_doa_noMC   ./ cnt_doa);

flr = @(x) max(x, 1e-3);
colorsM   = lines(nMethods);
markers_m = {'o-','s-','d-','^-'};

% --- (A) convergencia por SNR ---
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fig = figure('Name',sprintf('MultiDir SNR=%+d dB',SNR_dB),'Color','w','Position',[80 80 1250 520]);

    subplot(1,2,1); hold on; grid on;
    for im = 1:nMethods
        plot(1:maxIter, flr(squeeze(rmse_doa_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:),'LineWidth',1.6,'MarkerFaceColor',colorsM(im,:), ...
            'DisplayName',sprintf('%s (multi-dir)',methods{im}));
        yline(flr(rmse_doa_noComp(im,iSNR)), ':','Color',colorsM(im,:), ...
            'LineWidth',1.0,'HandleVisibility','off');
    end
    set(gca,'YScale','log'); xlabel('Iteracao'); ylabel('RMSE de DoA (graus)');
    title(sprintf('RMSE de DoA vs iteracao  (SNR=%+d dB, P=%d)',SNR_dB,P_default));
    legend('Location','best'); xlim([1 maxIter]);
    subtitle('linha pontilhada (...) = RMSE sem compensacao');

    subplot(1,2,2); hold on; grid on;
    for im = 1:nMethods
        plot(1:maxIter, flr(squeeze(mean_frob_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:),'LineWidth',1.6,'MarkerFaceColor',colorsM(im,:), ...
            'DisplayName',methods{im});
    end
    set(gca,'YScale','log'); xlabel('Iteracao');
    ylabel('media de ||C_{hat}-C_{true}||_F / ||C_{true}||_F');
    title('Erro da matriz de acoplamento vs iteracao');
    legend('Location','best'); xlim([1 maxIter]);

    exportgraphics(fig, fullfile(outDir, sprintf('multidir_iter_SNR_%+03d.png',SNR_dB)),'Resolution',180);
end

% --- (C) resumo vs SNR ---
fig = figure('Color','w','Position',[100 100 920 560]); hold on; grid on;
for im = 1:nMethods
    plot(range_SNR_dB, flr(rmse_doa_final(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',2.0,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8, ...
        'DisplayName',sprintf('%s multi-dir',methods{im}));
    plot(range_SNR_dB, flr(rmse_doa_noComp(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.2,'LineStyle',':','MarkerFaceColor','none','MarkerSize',8, ...
        'DisplayName',sprintf('%s sem comp.',methods{im}));
end
set(gca,'YScale','log'); xlabel('SNR (dB)'); ylabel('RMSE de DoA (graus)');
xticks(range_SNR_dB); xlim([min(range_SNR_dB)-1, max(range_SNR_dB)+1]);
title(sprintf(['RMSE de DoA: multi-dir (solido) vs sem comp (pontilhado)\n' ...
               'P=%d direcoes, %d rodadas, raio=%.2f\\lambda'], P_default, n_runs, r/lambda));
legend('Location','eastoutside');
exportgraphics(fig, fullfile(outDir,'resumo_RMSE_doa_vs_snr.png'),'Resolution',180);

%% =====================================================================
%   PARTE B: varredura do numero de direcoes P  (SNR = snr_Psweep)
%% =====================================================================
fprintf('\n===== [B] Varredura de P em SNR=%+d dB =====\n', snr_Psweep);
nP = numel(range_P);
finalRMSE_vsP = zeros(nMethods, nP);
finalFrob_vsP = zeros(nMethods, nP);
for iP = 1:nP
    P = range_P(iP);
    sumsq_doa_P = zeros(nMethods,1);  cnt_doa_P = 0;
    sum_frob_P  = zeros(nMethods,1);  cnt_run_P = 0;
    for run = 1:n_runs
        [B, Xc_cell, q_cell, phitrue] = gen_shots(P, grid_step, min_sep_deg, ...
            M, r, lambda, theta_sig_deg, snr_Psweep, N, fs, Rs, sps, alpha, span, fd, C_true);
        for im = 1:nMethods
            [frob_it, ssq_it] = multidir_selfcal_one(methods{im}, B, Xc_cell, q_cell, ...
                phitrue, M, r, lambda, theta_sig_deg, beta_uca, A_dict, phi_grid_deg, ...
                maxIter, C_true, normCtrue, inner_iter, inner_tol);
            sumsq_doa_P(im) = sumsq_doa_P(im) + ssq_it(end);
            sum_frob_P(im)  = sum_frob_P(im)  + frob_it(end);
        end
        cnt_doa_P = cnt_doa_P + P;
        cnt_run_P = cnt_run_P + 1;
    end
    finalRMSE_vsP(:,iP) = sqrt(sumsq_doa_P / cnt_doa_P);
    finalFrob_vsP(:,iP) = sum_frob_P / cnt_run_P;
    fprintf('  P=%d feito (%.0fs)\n', P, toc(t0));
end

fig = figure('Color','w','Position',[100 100 1250 520]);
subplot(1,2,1); hold on; grid on;
for im = 1:nMethods
    plot(range_P, flr(finalRMSE_vsP(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',methods{im});
end
set(gca,'YScale','log'); xlabel('nº de direcoes P'); ylabel('RMSE de DoA final (graus)');
xticks(range_P); title(sprintf('RMSE de DoA final vs P  (SNR=%+d dB)', snr_Psweep));
legend('Location','best');

subplot(1,2,2); hold on; grid on;
for im = 1:nMethods
    plot(range_P, flr(finalFrob_vsP(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',methods{im});
end
set(gca,'YScale','log'); xlabel('nº de direcoes P');
ylabel('Frobenius final ||C_{hat}-C_{true}||_F / ||C_{true}||_F');
xticks(range_P); title('Erro final de C vs P');
legend('Location','best');
sgtitle(sprintf('Ganho de identificabilidade com mais direcoes  (SNR=%+d dB, %d rodadas)', ...
    snr_Psweep, n_runs), 'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'varredura_numero_direcoes_P.png'),'Resolution',180);

%% ---- Resumo numerico ----
fprintf('\n========== RESUMO multi-dir (RMSE DoA, graus), P=%d ==========\n', P_default);
for iSNR = 1:nSNR
    fprintf('SNR = %+3d dB:\n', range_SNR_dB(iSNR));
    for im = 1:nMethods
        fprintf('   %-6s  multiDir=%8.3f   semComp=%8.3f   semAcop=%8.3f\n', ...
            methods{im}, rmse_doa_final(im,iSNR), rmse_doa_noComp(im,iSNR), rmse_doa_noMC(im,iSNR));
    end
end
fprintf('\nFiguras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function [B, Xc_cell, q_cell, phitrue] = gen_shots(P, grid_step, min_sep, ...
    M, r, lambda, theta_deg, SNR_dB, N, fs, Rs, sps, alpha, span, fd, C_true)
% GEN_SHOTS  Sorteia P direcoes distintas (no grid, separacao minima) e gera
%   os P shots NAO-SIMULTANEOS (um sinal por vez), retornando assinaturas B,
%   dados acoplados Xc_cell, struct q_cell (q + X_ideal) e direcoes verdadeiras.
    phitrue = zeros(1, P);
    phitrue(1) = round((-180 + 360*rand)/grid_step)*grid_step;
    for p = 2:P
        ok = false;
        while ~ok
            cand = round((-180 + 360*rand)/grid_step)*grid_step;
            if all(abs(wrapTo180(cand - phitrue(1:p-1))) >= min_sep)
                phitrue(p) = cand; ok = true;
            end
        end
    end
    B = zeros(M, P); Xc_cell = cell(1,P); q_cell = cell(1,P);
    for p = 1:P
        [~, q, ~, Xsig, ~, Xn] = utils.simulate_fsk_data_uca(M, r, lambda, ...
            phitrue(p), phitrue(p)+90, theta_deg, theta_deg, ...
            SNR_dB, -100, N, fs, Rs, sps, alpha, span, fd);
        q = q(:);
        Xc = C_true*Xsig + Xn;          % acoplado (ruido no receptor)
        Xc_cell{p} = Xc;
        s.q = q; s.X_ideal = Xsig + Xn; % ideal (sem acoplamento)
        q_cell{p} = s;
        B(:,p) = Xc*conj(q)/(q'*q);     % assinatura acoplada (fixa)
    end
end

function [frob_it, ssq_it] = multidir_selfcal_one(method, B, Xc_cell, q_cell, ...
    phitrue, M, r, lambda, theta_deg, beta_uca, A_dict, phi_grid_deg, ...
    maxIter, C_true, normCtrue, inner_iter, inner_tol)
% MULTIDIR_SELFCAL_ONE  Laco alternante multi-dir para UM metodo de DoA.
%   init = KW; C conjunto via estimate_C_multidir; sem damping/parada.
%   Retorna Frobenius por iteracao e SOMA de (erro DoA)^2 por iteracao (sobre P).
    P = numel(Xc_cell);

    % ----- Init mode = KW (direcoes iniciais via KW cru por shot) -----
    A0 = zeros(M, P);
    for p = 1:P
        phi0 = doa_estimate(Xc_cell{p}, 'KW', q_cell{p}.q, r, lambda, beta_uca, A_dict, phi_grid_deg);
        A0(:,p) = utils.steering_vec_uca(M, r, lambda, theta_deg, phi0);
    end
    C = estimate_C_multidir(B, A0, M, inner_iter, inner_tol);

    frob_it = zeros(maxIter,1);
    ssq_it  = zeros(maxIter,1);
    for it = 1:maxIter
        D = inv(C);
        A = zeros(M, P);
        ssq = 0;
        for p = 1:P
            Y   = D * Xc_cell{p};
            phi = doa_estimate(Y, method, q_cell{p}.q, r, lambda, beta_uca, A_dict, phi_grid_deg);
            A(:,p) = utils.steering_vec_uca(M, r, lambda, theta_deg, phi);
            ssq = ssq + wrapTo180(phi - phitrue(p))^2;
        end
        C = estimate_C_multidir(B, A, M, inner_iter, inner_tol);
        frob_it(it) = norm(C - C_true,'fro')/normCtrue;
        ssq_it(it)  = ssq;
    end
end

function phi_hat = doa_estimate(X, method, q, r, lambda, beta_uca, A_dict, phi_grid_deg)
% DOA_ESTIMATE  Estima o azimute (graus) por um dos 4 metodos.
    M = size(X,1);
    switch upper(method)
        case 'KW'
            [~, phi] = doa_kw_uca(X, q(:).', r, lambda, beta_uca);
            phi_hat  = phi(1);
        case 'DAS'
            R      = (X*X')/size(X,2);
            scores = real(sum(conj(A_dict).*(R*A_dict), 1));
            [~,ip] = max(scores); phi_hat = phi_grid_deg(ip);
        case 'CAPON'
            R      = (X*X')/size(X,2);
            R      = R + 1e-6*trace(R)/M*eye(M);
            Rinv   = R\eye(M);
            den    = real(sum(conj(A_dict).*(Rinv*A_dict), 1));
            scores = 1./max(den, eps); [~,ip] = max(scores); phi_hat = phi_grid_deg(ip);
        case 'MUSIC'
            R = (X*X')/size(X,2); R = (R+R')/2;
            [V, Dg]  = eig(R); [~, idx] = sort(real(diag(Dg)), 'descend'); V = V(:, idx);
            En   = V(:, 2:end);           % 1 fonte (sem interferente)
            proj = En'*A_dict;
            den  = real(sum(conj(proj).*proj, 1));
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
