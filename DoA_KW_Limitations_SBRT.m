% =========================================================================
%  DoA_KW_Limitations_SBRT.m
%
%  Limitações práticas do estimador de DOA com forma de onda conhecida (KW).
%  Acompanha o artigo SBrT 2026 "Practical Limitations of the Known Waveform
%  DOA Estimation Algorithm" e estende a simulação usada no EUSIPCO 2026.
%
%  Três estudos independentes, comparando KW com DAS, Capon e MUSIC:
%    (1) RMSE vs ISR para diferentes modulações do interferidor
%    (2) RMSE vs erro de sincronização (inteiro e fracionário)
%    (3) RMSE vs tamanho da sequência de treinamento (snapshots K)
%
%  Cada ponto é uma média Monte Carlo sobre `nEnsembles` realizações nas
%  quais bits, símbolos, modulantes, fase inicial e ruído são regerados.
%
%  Dependências:
%    - utils.m       (com simulate_data_uca_v2, gen_interferer, etc.)
%    - doa_kw_uca.m  (mesmo do EUSIPCO)
%    - MATLAB Antenna Toolbox (para compute_Ctx_for_R)
%
%  Autor: Nicolas S. M. M. de Oliveira  |  IME/SE-3
% =========================================================================

clear; clc; close all;

%% ===================== Parâmetros gerais ===============================
M       = 8;             % nº de elementos do UCA
fc      = 500e6;         % Hz
c       = 3e8;
lambda  = c/fc;
radius  = 0.25*lambda;   % raio compacto (mesmo do EUSIPCO)

theta_sig_deg = 90;      % plano XY
theta_int_deg = 90;

fs      = 288000;        % Hz - taxa de amostragem
N       = 21000;         % comprimento total disponível
N_DOA   = 8000;          % máximo de snapshots usados para DOA

% Parâmetros da forma de onda do SOI (2-FSK do EUSIPCO)
Rs      = 9600;
sps     = 30;
alpha   = 0.3;
span    = 8;
fd      = 4.8e3;

% Ângulos
phi_sig_deg_fixed = 30;
phi_int_deg_fixed = -45;
randomize_angles  = true;   % sorteia ângulos a cada ensemble

% Acoplamento mútuo
Z0 = 10;
apply_coupling = true;
fprintf('Calculando matriz de acoplamento para r=%.3f lambda...\n', radius/lambda);
Ctx_full = utils.compute_Ctx_for_R(fc, M, radius, Z0);
CM_id    = eye(M);
if apply_coupling
    CM = Ctx_full;
else
    CM = CM_id;
end

% Tipos de interferidor estudados
interferer_types = ["FM","FSK2","QAM64","NOISE"];

% Monte Carlo
nEnsembles = 100;     % nº de realizações por ponto (sobe para ~50 no run final)

% rng(42);  % descomente para reprodutibilidade

% Pasta de saída
outDir = fullfile(pwd, 'Graficos_SBRT_KW_Limitations');
if ~exist(outDir, 'dir'); mkdir(outDir); end

%% ===================== Métodos e estilo ================================
method_names = ["KW", "DAS", "Capon", "MUSIC"];
nMethods     = numel(method_names);

% Mesma paleta do paper EUSIPCO para consistência visual
method_markers = ["o-","x-","s-","d-"];
method_colors  = ["red", "green", "blue", "black"];

%% ===================== Pré-cálculo do scan angular =====================
% Para os métodos clássicos (DAS, Capon, MUSIC). Calculado UMA VEZ no script.
phi_scan = -180:0.05:180;          % graus (0.05° => 7201 pontos)
Ngrid    = numel(phi_scan);
A_scan   = zeros(M, Ngrid);
for ig = 1:Ngrid
    A_scan(:,ig) = utils.steering_vec_uca(M, radius, lambda, ...
                                          theta_sig_deg, phi_scan(ig));
end
fprintf('Scan grid: %d pontos (resolução %.3f°)\n', Ngrid, phi_scan(2)-phi_scan(1));

%% ===================== Wrapper único ===================================
% Roda UMA estimação por cada um dos 4 métodos para um cenário completo.
% Retorna um vetor 4x1 com [phi_KW; phi_DAS; phi_Capon; phi_MUSIC] em graus.
run_all_once = @(SNR_dB, ISR_dB, interferer_type, sync_offset, K, ...
                 phi_sig_deg, phi_int_deg) ...
    do_one_estimation_all(M, radius, lambda, CM, ...
                          phi_sig_deg, phi_int_deg, ...
                          theta_sig_deg, theta_int_deg, ...
                          SNR_dB, ISR_dB, N, fs, ...
                          Rs, sps, alpha, span, fd, ...
                          interferer_type, sync_offset, K, ...
                          A_scan, phi_scan);

%% =======================================================================
%  ESTUDO 1: RMSE vs ISR para diferentes modulações do interferidor
% =======================================================================
% fprintf('\n========================================================\n');
% fprintf('Estudo 1: RMSE vs ISR (KW, DAS, Capon, MUSIC)\n');
% fprintf('========================================================\n');
% 
% SNR_dB_S1   = 10;
% K_S1        = 2000;
% sync_S1     = 0;
% range_ISR_S1 = -18:3:6;
% 
% nMod = numel(interferer_types);
% nISR = numel(range_ISR_S1);
% RMSE_S1 = zeros(nMethods, nMod, nISR);
% 
% t0 = tic;
% for iMod = 1:nMod
%     for iISR = 1:nISR
%         errs = zeros(nMethods, nEnsembles);
%         for e = 1:nEnsembles
%             if randomize_angles
%                 phi_s = -180 + 360*rand;
%                 phi_i = -180 + 360*rand;
%             else
%                 phi_s = phi_sig_deg_fixed;
%                 phi_i = phi_int_deg_fixed;
%             end
%             phi_hats = run_all_once(SNR_dB_S1, range_ISR_S1(iISR), ...
%                                     interferer_types(iMod), sync_S1, K_S1, ...
%                                     phi_s, phi_i);
%             errs(:,e) = abs(mod(phi_hats - phi_s + 180, 360) - 180);
%         end
%         RMSE_S1(:,iMod,iISR) = sqrt(mean(errs.^2, 2));
%         fprintf('  [Mod=%-6s | ISR=%+3d dB] RMSE: KW=%6.2f  DAS=%6.2f  Capon=%6.2f  MUSIC=%6.2f  (%.1fs)\n', ...
%                 interferer_types(iMod), range_ISR_S1(iISR), ...
%                 RMSE_S1(1,iMod,iISR), RMSE_S1(2,iMod,iISR), ...
%                 RMSE_S1(3,iMod,iISR), RMSE_S1(4,iMod,iISR), toc(t0));
%     end
% end
% 
% % Gráficos S1: uma figura por modulação, 4 curvas (uma por método)
% for iMod = 1:nMod
%     figname = sprintf('S1 RMSE vs ISR | Interf = %s', interferer_types(iMod));
%     fig = figure('Name', figname, 'NumberTitle','off', ...
%                  'Position',[100 100 800 500]);
%     hold on; grid on; box on;
%     for m = 1:nMethods
%         plot(range_ISR_S1, squeeze(RMSE_S1(m,iMod,:)), ...
%              method_markers(m), 'Color', method_colors(m), ...
%              'LineWidth', 1.6, 'MarkerSize', 7);
%     end
%     set(gca,'YScale','log');
%     xlabel('ISR (dB)','Interpreter','latex');
%     ylabel('RMSE (degrees)','Interpreter','latex');
%     ylim([0 110]);
%     title(sprintf('RMSE vs ISR | Interf = %s (SNR=%d dB, K=%d)', ...
%                   interferer_types(iMod), SNR_dB_S1, K_S1));
%     legend(method_names,'Location','best','Interpreter','latex');
% 
%     base = sprintf('S1_RMSE_vs_ISR_Mod_%s', interferer_types(iMod));
%     exportgraphics(fig, fullfile(outDir,[char(base) '.png']),'Resolution',300);
%     try
%         matlab2tikz(fullfile(outDir,[char(base) '.tex']), ...
%                     'width','\figurewidth','height','\figureheight');
%     catch ME
%         warning('matlab2tikz indisponível: %s', ME.message);
%     end
% end
% 
% % --- Gráfico resumo (KW apenas, todas as modulações) ---
% fig = figure('Name','S1 KW: RMSE vs ISR (todas as modulações)', ...
%              'NumberTitle','off','Position',[100 100 800 500]);
% hold on; grid on; box on;
% summary_markers = ["o-","s-","d-","^-","v-","*-","x-","+-"];
% summary_colors  = lines(nMod);
% for iMod = 1:nMod
%     plot(range_ISR_S1, squeeze(RMSE_S1(1,iMod,:)), ...
%          summary_markers(min(iMod,numel(summary_markers))), ...
%          'Color', summary_colors(iMod,:), ...
%          'LineWidth', 1.6, 'MarkerSize', 7);
% end
% set(gca,'YScale','log');
% xlabel('ISR (dB)','Interpreter','latex');
% ylabel('RMSE (degrees)','Interpreter','latex');
% ylim([0 110]);
% title(sprintf('KW: RMSE vs ISR (SNR=%d dB, K=%d, sync=0)', SNR_dB_S1, K_S1));
% legend(interferer_types,'Location','best','Interpreter','latex');
% exportgraphics(fig, fullfile(outDir,'S1_KW_RMSE_vs_ISR_by_Modulation.png'),'Resolution',300);
% try
%     matlab2tikz(fullfile(outDir,'S1_KW_RMSE_vs_ISR_by_Modulation.tex'), ...
%                 'width','\figurewidth','height','\figureheight');
% catch ME
%     warning('matlab2tikz indisponível: %s', ME.message);
% end
% 
% save(fullfile(outDir,'results_S1.mat'), 'RMSE_S1','range_ISR_S1', ...
%      'interferer_types','SNR_dB_S1','K_S1','method_names');
% 
% return;

%% =======================================================================
%  ESTUDO 2: RMSE vs erro de sincronização
% =======================================================================
% fprintf('\n========================================================\n');
% fprintf('Estudo 2: RMSE vs erro de sincronização (KW, DAS, Capon, MUSIC)\n');
% fprintf('========================================================\n');
% 
% % Erro de sincronização expresso como FRAÇÃO do período de símbolo T_sym.
% % T_sym = sps amostras (=30 amostras a fs=288 kHz). A autocorrelação de uma
% % 2-FSK colapsa em ~1 T_sym, então faz sentido amostrar densamente
% % sub-T_sym e estender até alguns T_sym para capturar o "joelho" da curva.
% % Internamente, delayseq recebe o valor em amostras (suporta fracionário).
% range_sync_Tsym = [-30 -20 -10 -4 -3 -2 -1 -0.5 ...
%                     0 ...
%                     0.5 1 2 3 4 10 20 30];
% range_sync = range_sync_Tsym * sps;   % converte para amostras
% 
% % sync_scenarios = struct( ...
% %     'SNR_dB', { 0,  6, 12, 12}, ...
% %     'ISR_dB', {-20, -20, -20, -3});
% sync_scenarios = struct( ...
%     'SNR_dB', { 12}, ...
%     'ISR_dB', {-20});
% 
% nScen = numel(sync_scenarios);
% K_S2  = 2000;
% mod_S2 = "FM";
% 
% nSync = numel(range_sync);
% RMSE_S2 = zeros(nMethods, nScen, nSync);
% 
% t0 = tic;
% for iScen = 1:nScen
%     SNR_dB = sync_scenarios(iScen).SNR_dB;
%     ISR_dB = sync_scenarios(iScen).ISR_dB;
%     for iSync = 1:nSync
%         errs = zeros(nMethods, nEnsembles);
%         for e = 1:nEnsembles
%             if randomize_angles
%                 phi_s = -180 + 360*rand;
%                 phi_i = -180 + 360*rand;
%             else
%                 phi_s = phi_sig_deg_fixed;
%                 phi_i = phi_int_deg_fixed;
%             end
%             phi_hats = run_all_once(SNR_dB, ISR_dB, mod_S2, ...
%                                     range_sync(iSync), K_S2, ...
%                                     phi_s, phi_i);
%             errs(:,e) = abs(mod(phi_hats - phi_s + 180, 360) - 180);
%         end
%         RMSE_S2(:,iScen,iSync) = sqrt(mean(errs.^2, 2));
%         fprintf('  [SNR=%+3d ISR=%+3d sync=%+5.2f T_sym (%+6.2f samp)] RMSE: KW=%6.2f  DAS=%6.2f  Capon=%6.2f  MUSIC=%6.2f  (%.1fs)\n', ...
%                 SNR_dB, ISR_dB, range_sync_Tsym(iSync), range_sync(iSync), ...
%                 RMSE_S2(1,iScen,iSync), RMSE_S2(2,iScen,iSync), ...
%                 RMSE_S2(3,iScen,iSync), RMSE_S2(4,iScen,iSync), toc(t0));
%     end
% end
% 
% % Gráficos S2: uma figura por cenário SNR/ISR, 4 curvas (uma por método)
% for iScen = 1:nScen
%     SNR_dB = sync_scenarios(iScen).SNR_dB;
%     ISR_dB = sync_scenarios(iScen).ISR_dB;
%     figname = sprintf('S2 RMSE vs Sync | SNR=%+d ISR=%+d', SNR_dB, ISR_dB);
%     fig = figure('Name', figname, 'NumberTitle','off', ...
%                  'Position',[100 100 800 500]);
%     hold on; grid on; box on;
%     for m = 1:nMethods
%         plot(range_sync_Tsym, squeeze(RMSE_S2(m,iScen,:)), ...
%              method_markers(m), 'Color', method_colors(m), ...
%              'LineWidth', 1.6, 'MarkerSize', 6);
%     end
%     set(gca,'YScale','log');
%     xlabel('Erro de sincronização ($\Delta / T_{\mathrm{sym}}$)','Interpreter','latex');
%     ylabel('RMSE (degrees)','Interpreter','latex');
%     title(sprintf('RMSE vs sync error | SNR=%+d, ISR=%+d (Interf=%s, K=%d)', ...
%                   SNR_dB, ISR_dB, mod_S2, K_S2));
%     legend(method_names,'Location','best','Interpreter','latex');
% 
%     base = sprintf('S2_RMSE_vs_Sync_SNR_%+ddB_ISR_%+ddB', SNR_dB, ISR_dB);
%     base = strrep(base,'+','p');  base = strrep(base,'-','m');
%     exportgraphics(fig, fullfile(outDir,[base '.png']),'Resolution',300);
%     try
%         matlab2tikz(fullfile(outDir,[base '.tex']), ...
%                     'width','\figurewidth','height','\figureheight');
%     catch ME
%         warning('matlab2tikz indisponível: %s', ME.message);
%     end
% end
% 
% % --- Gráfico resumo (KW apenas, todos os cenários SNR/ISR) ---
% fig = figure('Name','S2 KW: RMSE vs sync (todos os cenários)', ...
%              'NumberTitle','off','Position',[100 100 800 500]);
% hold on; grid on; box on;
% summary_markers = ["o-","s-","d-","^-","v-","*-","x-","+-"];
% summary_colors  = lines(nScen);
% legend_strings = strings(nScen,1);
% for iScen = 1:nScen
%     plot(range_sync_Tsym, squeeze(RMSE_S2(1,iScen,:)), ...
%          summary_markers(min(iScen,numel(summary_markers))), ...
%          'Color', summary_colors(iScen,:), ...
%          'LineWidth', 1.6, 'MarkerSize', 6);
%     legend_strings(iScen) = sprintf('SNR=%+d, ISR=%+d', ...
%         sync_scenarios(iScen).SNR_dB, sync_scenarios(iScen).ISR_dB);
% end
% set(gca,'YScale','log');
% xlabel('Erro de sincronização ($\Delta / T_{\mathrm{sym}}$)','Interpreter','latex');
% ylabel('RMSE (degrees)','Interpreter','latex');
% title(sprintf('KW: RMSE vs sync error (Interf=%s, K=%d)', mod_S2, K_S2));
% legend(legend_strings,'Location','best','Interpreter','latex');
% exportgraphics(fig, fullfile(outDir,'S2_KW_RMSE_vs_SyncOffset.png'),'Resolution',300);
% try
%     matlab2tikz(fullfile(outDir,'S2_KW_RMSE_vs_SyncOffset.tex'), ...
%                 'width','\figurewidth','height','\figureheight');
% catch ME
%     warning('matlab2tikz indisponível: %s', ME.message);
% end
% 
% save(fullfile(outDir,'results_S2.mat'), 'RMSE_S2','range_sync','range_sync_Tsym', ...
%      'sync_scenarios','K_S2','mod_S2','method_names');
% 
% return;

%% =======================================================================
%  ESTUDO 3: RMSE vs tamanho da sequência de treinamento (snapshots)
% =======================================================================
fprintf('\n========================================================\n');
fprintf('Estudo 3: RMSE vs número de snapshots K (KW, DAS, Capon, MUSIC)\n');
fprintf('========================================================\n');

range_K  = [50 100 200 400 800 1500 2500 4000 6000 8000];
range_K  = range_K(range_K <= N_DOA);

K_scenarios = struct( ...
    'SNR_dB', {-3,  0,  6,  6}, ...
    'ISR_dB', {-6,  0, -6,  0});
nKScen = numel(K_scenarios);
sync_S3 = 0;
mod_S3  = "FM";

nK = numel(range_K);
RMSE_S3 = zeros(nMethods, nKScen, nK);

t0 = tic;
for iScen = 1:nKScen
    SNR_dB = K_scenarios(iScen).SNR_dB;
    ISR_dB = K_scenarios(iScen).ISR_dB;
    for iK = 1:nK
        K = range_K(iK);
        errs = zeros(nMethods, nEnsembles);
        for e = 1:nEnsembles
            if randomize_angles
                phi_s = -180 + 360*rand;
                phi_i = -180 + 360*rand;
            else
                phi_s = phi_sig_deg_fixed;
                phi_i = phi_int_deg_fixed;
            end
            phi_hats = run_all_once(SNR_dB, ISR_dB, mod_S3, sync_S3, K, ...
                                    phi_s, phi_i);
            errs(:,e) = abs(mod(phi_hats - phi_s + 180, 360) - 180);
        end
        RMSE_S3(:,iScen,iK) = sqrt(mean(errs.^2, 2));
        fprintf('  [SNR=%+3d ISR=%+3d K=%5d] RMSE: KW=%6.2f  DAS=%6.2f  Capon=%6.2f  MUSIC=%6.2f  (%.1fs)\n', ...
                SNR_dB, ISR_dB, K, ...
                RMSE_S3(1,iScen,iK), RMSE_S3(2,iScen,iK), ...
                RMSE_S3(3,iScen,iK), RMSE_S3(4,iScen,iK), toc(t0));
    end
end

% Gráficos S3: uma figura por cenário SNR/ISR, 4 curvas (uma por método)
for iScen = 1:nKScen
    SNR_dB = K_scenarios(iScen).SNR_dB;
    ISR_dB = K_scenarios(iScen).ISR_dB;
    figname = sprintf('S3 RMSE vs K | SNR=%+d ISR=%+d', SNR_dB, ISR_dB);
    fig = figure('Name', figname, 'NumberTitle','off', ...
                 'Position',[100 100 800 500]);
    hold on; grid on; box on;
    for m = 1:nMethods
        plot(range_K, squeeze(RMSE_S3(m,iScen,:)), ...
             method_markers(m), 'Color', method_colors(m), ...
             'LineWidth', 1.6, 'MarkerSize', 6);
    end
    set(gca,'YScale','log');
    set(gca,'XScale','log');
    xlabel('Número de snapshots K','Interpreter','latex');
    ylabel('RMSE (degrees)','Interpreter','latex');
    title(sprintf('RMSE vs K | SNR=%+d, ISR=%+d (Interf=%s, sync=0)', ...
                  SNR_dB, ISR_dB, mod_S3));
    legend(method_names,'Location','best','Interpreter','latex');

    base = sprintf('S3_RMSE_vs_K_SNR_%+ddB_ISR_%+ddB', SNR_dB, ISR_dB);
    base = strrep(base,'+','p');  base = strrep(base,'-','m');
    exportgraphics(fig, fullfile(outDir,[base '.png']),'Resolution',300);
    try
        matlab2tikz(fullfile(outDir,[base '.tex']), ...
                    'width','\figurewidth','height','\figureheight');
    catch ME
        warning('matlab2tikz indisponível: %s', ME.message);
    end
end

% --- Gráfico resumo (KW apenas, todos os cenários SNR/ISR) ---
fig = figure('Name','S3 KW: RMSE vs K (todos os cenários)', ...
             'NumberTitle','off','Position',[100 100 800 500]);
hold on; grid on; box on;
summary_markers = ["o-","s-","d-","^-","v-","*-","x-","+-"];
summary_colors  = lines(nKScen);
legend_strings = strings(nKScen,1);
for iScen = 1:nKScen
    plot(range_K, squeeze(RMSE_S3(1,iScen,:)), ...
         summary_markers(min(iScen,numel(summary_markers))), ...
         'Color', summary_colors(iScen,:), ...
         'LineWidth', 1.6, 'MarkerSize', 6);
    legend_strings(iScen) = sprintf('SNR=%+d, ISR=%+d', ...
        K_scenarios(iScen).SNR_dB, K_scenarios(iScen).ISR_dB);
end
set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel('Número de snapshots K','Interpreter','latex');
ylabel('RMSE (degrees)','Interpreter','latex');
title(sprintf('KW: RMSE vs K (Interf=%s, sync=0)', mod_S3));
legend(legend_strings,'Location','best','Interpreter','latex');
exportgraphics(fig, fullfile(outDir,'S3_KW_RMSE_vs_Snapshots.png'),'Resolution',300);
try
    matlab2tikz(fullfile(outDir,'S3_KW_RMSE_vs_Snapshots.tex'), ...
                'width','\figurewidth','height','\figureheight');
catch ME
    warning('matlab2tikz indisponível: %s', ME.message);
end

save(fullfile(outDir,'results_S3.mat'), 'RMSE_S3','range_K', ...
     'K_scenarios','mod_S3','sync_S3','method_names');

fprintf('\nConcluído. Resultados salvos em %s\n', outDir);


% =========================================================================
%  FUNÇÃO LOCAL: Estimação DOA por TODOS os métodos (KW, DAS, Capon, MUSIC)
% =========================================================================
function phi_hats = do_one_estimation_all(M, radius, lambda, CM, ...
                                          phi_sig_deg, phi_int_deg, ...
                                          theta_sig_deg, theta_int_deg, ...
                                          SNR_dB, ISR_dB, N, fs, ...
                                          Rs, sps, alpha, span, fd, ...
                                          interferer_type, sync_offset, K, ...
                                          A_scan, phi_scan)
%DO_ONE_ESTIMATION_ALL  Retorna [phi_KW; phi_DAS; phi_Capon; phi_MUSIC] em graus.
%
%   A_scan : matriz M x Ngrid com os steering vectors pré-calculados em
%            phi_scan (graus). Usada pelos métodos clássicos.

    % ----- Geração do sinal -----
    [~, ~, ~, Xsig, Xint, Xn, ~, ~, ~, q_local] = ...
        utils.simulate_data_uca_v2(M, radius, lambda, ...
                                   phi_sig_deg, phi_int_deg, ...
                                   theta_sig_deg, theta_int_deg, ...
                                   SNR_dB, ISR_dB, N, fs, ...
                                   Rs, sps, alpha, span, fd, ...
                                   interferer_type, sync_offset);

    % Modelo do EUSIPCO: X = CM*Xsig + CM*Xint + Xn
    X = CM*Xsig + CM*Xint + Xn;
    K = min(K, size(X,2));
    Xk = X(:,1:K);

    % ----- (1) KW -----
    beta = 2*pi*(0:M-1)'/M;
    [~, phi_KW] = doa_kw_uca(Xk, q_local(1:K), radius, lambda, beta);

    % ----- Matriz de covariância (compartilhada entre os clássicos) -----
    Rxx   = (Xk*Xk')/K;
    delta = 1e-3 * trace(Rxx)/M;
    Rxx_dl = Rxx + delta*eye(M);

    % Decomposição espectral (para MUSIC)
    [eigvec, eigval] = eig(Rxx_dl);
    [~, idx] = sort(diag(eigval), 'descend');
    E   = eigvec(:, idx);
    Ksrc = 1;                   % SOI + 1 interferidor
    En  = E(:, Ksrc+1:end);
    EnEnH = En*En';

    Rinv = inv(Rxx_dl);

    % ----- Scan vetorizado (DAS, Capon, MUSIC) -----
    %   DAS    : P = a' * Rxx   * a
    %   Capon  : P = 1/(a' * Rinv * a)
    %   MUSIC  : P = 1/(a' * En*En' * a)
    %
    % Calcula todos os a' * M * a de uma vez via
    %   sum(conj(A_scan) .* (M*A_scan), 1).
    P_DAS   = real(sum(conj(A_scan) .* (Rxx    * A_scan), 1));    % 1 x Ngrid
    P_Capon = 1 ./ max(real(sum(conj(A_scan) .* (Rinv  * A_scan), 1)), eps);
    P_MUSIC = 1 ./ max(real(sum(conj(A_scan) .* (EnEnH * A_scan), 1)), eps);

    [~, i_das]   = max(P_DAS);
    [~, i_capon] = max(P_Capon);
    [~, i_music] = max(P_MUSIC);

    phi_DAS    = phi_scan(i_das);
    phi_Capon  = phi_scan(i_capon);
    phi_MUSIC  = phi_scan(i_music);

    phi_hats = [phi_KW; phi_DAS; phi_Capon; phi_MUSIC];
end
