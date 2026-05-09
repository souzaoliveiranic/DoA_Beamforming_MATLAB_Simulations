% Antes esse código era phase-only, agora ele é delay+phase
clear; clc;
%rng(1);
close all;
%% Parâmetros do array
M      = 8;           % nº de elementos do ULA
fc     = 500e6; %2400e6         % Hz
c      = 3e8;
lambda = c/fc;
r = 0.25 * lambda;    % raio 1/4 λ

theta_sig_deg = 90;   % plano XY
theta_int_deg = 90;
sigma_erro_phi = 0; %10;

%nTrials  = 10; %360;   % nº de realizações
fs     = 288000;       % taxa de amostragem (Hz) para formar snapshots
N      = 2100;        % nº de amostras
N_DOA  = 2100;
t      = (0:N-1).' / fs;
N_fm   = 10000;      % nº de amostras do sinal fm

fm_msg = 1000;        % Hz - frequência da mensagem
delta_f   = 5000;     % desvio FM
phi0   = 2*pi*rand;   % fase inicial aleatória do sinal
epsilon = 0.1;
p = 1.1;
maxIter = 100;
eta = 1e-3;

Rs     = 9600;          % taxa de símbolos [sym/s]
sps    = 30;            % amostras por símbolo
alpha  = 0.3;           % roll-off do RRC
span   = 8;             % comprimento do RRC em símbolos (TX/RX)
fd     = 4.8e3;         % desvio de frequência (Δf) [Hz]

range_SNR_dB = 6; %-6:3:6;
range_ISR_dB = -1; %-6:1:-3;
range_snapshots = 2000:3000:N_DOA;
range_radius = [0.15]; %[0.25 0.2 0.15 0.1];
% Codigos de Coupling (renumerados):
%   0 = sem coupling (ideal)
%   1 = com coupling, SEM compensacao
%   2 = com coupling, comp. via modelo IDEAL  (D = inv(C_true))   "oracle"
%   3 = com coupling, comp. via modelo PERTURBADO (Z_t com erro aleatorio)
%   4 = com coupling, comp. via KW LS 1-direcao
%   5 = com coupling, comp. via KW LS varias-direcoes (P direcoes)
%   6 = com coupling, comp. via SELF-CAL com DAS
%   7 = com coupling, comp. via SELF-CAL com CAPON
%   8 = com coupling, comp. via SELF-CAL com MUSIC
%   9 = com coupling, comp. via SELF-CAL com KW
range_coupling = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];

% --- parametros da calibracao (Khan 2020): UMA medicao de uma direcao conhecida ---
phi_cal_deg   = 0;     % azimute da fonte de calibracao (Coupling=4)
theta_cal_deg = 90;    % elevacao (plano XY)
SNR_cal_dB    = 10;    % SNR alto (camara anecoica)
ISR_cal_dB    = -40;   % praticamente sem interferente

% --- Modelo perturbado (Coupling=3): incerteza relativa em Z_t ---
%       Simula erro de simulacao EM / variabilidade de fabricacao do PCB.
pert_level_rel = 0.15;     % 5% de incerteza relativa nos coeficientes
rng(2025, 'twister');      % reprodutibilidade

% --- KW varias direcoes (Coupling=5): P direcoes de calibracao ---
phi_cal_multi_deg = [0, 30, 60, 90];   % evita 22.5 (singular UCA-8)
P_multi = numel(phi_cal_multi_deg);

% --- Self-cal (Coupling=6..9): parametros do alternante ---
selfcal_max_iter = 12;
selfcal_grid_deg = -180:0.5:180;
selfcal_methods  = {'DAS', 'CAPON', 'MUSIC', 'KW'};   % ordem para Coupling 6..9
selfcal_damping  = 0.5;        % subrelaxacao: suaviza oscilacoes (1.0 = sem damping)
selfcal_init     = 'identity';       % 'identity' ou 'kw' (recomendado quando ha interferentes)
range_phi = [75 100]; %-180:18:180;
teste_phi = 75;

methodsDoa = ["KW","DAS","MPDR","MUSIC"];
methodsBeamforming = ["SEM BF","MPDR","DAS"];
markers = ["o-","x-","s-","d-","^-","v-","*-","+-"];
colors =  ["red", "green", "blue", "black", "magenta", "cyan", "yellow"];
line_style =   ["-", "--", ":","-."]; % linha continua sem acoplamento % linha tracejada com acoplamento

% Criar pasta 'graficos'
outDir = fullfile(pwd, 'new_graficos');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

nulo_no_interferidor = 0;

nMethodsDoa = numel(methodsDoa);
nMethodsBeamforming = numel(methodsBeamforming);
nSNR = numel(range_SNR_dB);
nISR = numel(range_ISR_dB);
nRadius  = numel(range_radius);
nSnapshots  = numel(range_snapshots);
nPhi  = numel(range_phi);
nCoupling  = numel(range_coupling);

TotalSim = nSNR*nISR*nRadius*nCoupling*nSnapshots*nPhi*nPhi;
iTotal = 0;
% Armazenamento
RMSE = zeros(nMethodsDoa, nSNR, nISR, nRadius, nCoupling, nSnapshots, nPhi*nPhi);
BER = zeros(4, nSNR, nISR, nRadius, nMethodsDoa, nPhi*nPhi);
EVM = zeros(4, nSNR, nISR, nRadius, nMethodsDoa, nPhi*nPhi);

%% ======= Coupling Matrix Estimation =======

Z0 = 50; % Impedância de referência

% Vetor com os valores únicos de acoplamento (da tabela) para 2.4GHz e
% 0.5*lambada
% Ordem: distância 1 até 7
Zt_vals = [ ...
   -0.08 - 1j*11.77;   % Z1
    7.19 + 1j*2.12;    % Z2
   -0.06 + 1j*5.96;    % Z3
   -2.50 + 1j*5.00;    % Z4
    0.00 + 1j*5.97;    % Z5
    7.26 + 1j*1.98;    % Z6
   -0.37 - 1j*11.85    % Z7
];

% Inicializa matriz Z
Z = eye(M);

% Preenche matriz
for i = 1:M
    for j = 1:M
        if i ~= j
            % Distância circular no UCA
            d = mod(abs(i - j), M);
            if d == 0
                d = M;
            end
            
            % Garante menor distância circular
            d = min(d, M - d);
            
            % MATLAB index começa em 1
            Z(i,j) = - Zt_vals(d) / Z0;
        end
    end
end

% Exibe resultado
disp(Z)

% C_true = matriz de acoplamento verdadeira (igual a Z aqui).
% (no codigo original havia um loop que acabava deixando C_true = eye(M); isso
%  era um bug porque a metrica de erro de Frobenius ficava sem sentido).
C_true = Z;

% Coeficientes unicos da estrutura circulante (para diagnostico e plot)
% Para M=8: c1, c2, c3 aparecem em 2 posicoes; c4 (distancia M/2) aparece em 1.
c_true_vec = zeros(M/2, 1);
for k = 1:M/2
    c_true_vec(k) = C_true(1, 1 + k);   % primeira linha, posicao 1+k
end


Coupling_matrices = zeros(M, M, nRadius); 
% Matrizes a ser usada nas simulações
for iRadius = 1:nRadius
    radius = range_radius(iRadius)*lambda;
    % Ctx = Z;
    Ctx = compute_Ctx_for_R(fc, M, radius, Z0);
    Coupling_matrices(:,:,iRadius) = Ctx; %inv(Ctx);
end

%% ============================================================================
%   CALIBRACAO (Khan 2020): UMA medicao numa direcao conhecida basta para
%   estimar a matriz de decoupling, porque o acoplamento depende SO da
%   geometria do array.
%
%   Modelo: b_hat_cal = alpha * C * a(theta_cal, phi_cal) + ruido
%   Estimacao: LS linear explorando a estrutura circulante simetrica do UCA.
%   Cada raio tem sua propria D_hat (na simulacao atual, range_radius=[0.5]).
%% ============================================================================

D_model_a    = cell(nRadius, 1);   % (Coupling=2) D oraculo: inv(C_true)
D_model_b    = cell(nRadius, 1);   % (Coupling=3) D modelo perturbado
D_kw_1dir    = cell(nRadius, 1);   % (Coupling=4) D via KW LS 1-direcao
D_kw_multi   = cell(nRadius, 1);   % (Coupling=5) D via KW LS varias-direcoes
C_kw_1dir    = cell(nRadius, 1);
C_model_a    = cell(nRadius, 1);
C_model_b    = cell(nRadius, 1);
C_kw_multi   = cell(nRadius, 1);
err_C_kw1    = zeros(nRadius, 1);  % KW 1-dir vs C_true
err_Cb_F     = zeros(nRadius, 1);  % modelo perturbado
err_C_kwmd   = zeros(nRadius, 1);  % KW multi-dir vs C_true
c_hat_vec    = cell(nRadius, 1);
c_pert_vec   = cell(nRadius, 1);
c_md_vec     = cell(nRadius, 1);

for iRadius = 1:nRadius
    radius_cal = range_radius(iRadius)*lambda;
    fprintf('\n=== Calibracao (raio = %.2f lambda) ===\n', range_radius(iRadius));

    % --- 1. Simula o cenario de calibracao (alta SNR, sem interferente) ---
    [~, q_cal, ~, Xsig_cal, Xint_cal, Xn_cal, ~, ~, ~, ~, ~, ~, ~] = ...
        utils.simulate_fsk_data_uca(M, radius_cal, lambda, ...
            phi_cal_deg, phi_cal_deg + 90, ...    % phi_int arbitrario, ja que ISR e' muito baixa
            theta_cal_deg, theta_cal_deg, ...
            SNR_cal_dB, ISR_cal_dB, N, fs, Rs, sps, alpha, span, fd);

    % Aplica acoplamento verdadeiro (mesmo modelo do loop principal)
    X_cal = Coupling_matrices(:,:,iRadius) * Xsig_cal + Coupling_matrices(:,:,iRadius) * Xint_cal;
    X_cal = X_cal ./ sqrt(mean(abs(X_cal).^2));
    X_cal = X_cal + Xn_cal;

    % --- 2. Estimativa da assinatura espacial via formas de onda conhecidas ---
    %  b_hat = X * q* / (q'*q)   (mesma logica do KW-DoA)
    b_hat_cal = X_cal * conj(q_cal(:)) / (q_cal(:)' * q_cal(:));

    % --- 3. Vetor de direcao ideal na direcao de calibracao ---
    a_cal = utils.steering_vec_uca(M, radius_cal, lambda, theta_cal_deg, phi_cal_deg);

    % --- 4. (Coupling=4) KW LS 1-direcao ---------------------------------
    [C_kw1_iR, c_hat_iR, alpha_hat_iR, residual_iR] = ...
        estimate_C_circulant_uca(b_hat_cal, a_cal, M);

    C_kw_1dir{iRadius} = C_kw1_iR;
    D_kw_1dir{iRadius} = inv(C_kw1_iR);
    err_C_kw1(iRadius) = norm(C_kw1_iR - C_true, 'fro') / norm(C_true, 'fro');
    c_hat_vec{iRadius} = c_hat_iR;

    % --- (Coupling=2) Modelo IDEAL: D = inv(C_true) (oracle) -------------
    C_model_a{iRadius} = C_true;
    D_model_a{iRadius} = inv(C_true);

    % --- (Coupling=3) Modelo PERTURBADO: Z_t com erro relativo gaussiano -
    Zt_pert = Zt_vals .* ( 1 + pert_level_rel * ...
              ( randn(size(Zt_vals)) + 1j*randn(size(Zt_vals)) ) / sqrt(2) );
    c_pert  = -Zt_pert / Z0;
    first_row_pert = [1, c_pert(1), c_pert(2), c_pert(3), ...
                         c_pert(4), ...
                         c_pert(3), c_pert(2), c_pert(1)];
    C_pert  = zeros(M);
    for ii_p = 1:M
        C_pert(ii_p, :) = circshift(first_row_pert, [0, ii_p-1]);
    end
    C_model_b{iRadius}  = C_pert;
    D_model_b{iRadius}  = inv(C_pert);
    err_Cb_F(iRadius)   = norm(C_pert - C_true, 'fro') / norm(C_true, 'fro');
    c_pert_vec{iRadius} = c_pert;

    % --- (Coupling=5) KW LS varias-direcoes ------------------------------
    B_multi = zeros(M, P_multi);
    A_multi = zeros(M, P_multi);
    for pp = 1:P_multi
        [~, q_p, ~, Xs_p, Xi_p, Xn_p, ~,~,~,~,~,~,~] = ...
            utils.simulate_fsk_data_uca(M, radius_cal, lambda, ...
                phi_cal_multi_deg(pp), phi_cal_multi_deg(pp) + 90, ...
                theta_cal_deg, theta_cal_deg, ...
                SNR_cal_dB, ISR_cal_dB, N, fs, Rs, sps, alpha, span, fd);
        X_p = Coupling_matrices(:,:,iRadius) * Xs_p + ...
              Coupling_matrices(:,:,iRadius) * Xi_p;
        X_p = X_p ./ sqrt(mean(abs(X_p).^2));
        X_p = X_p + Xn_p;
        B_multi(:, pp) = X_p * conj(q_p(:)) / (q_p(:)' * q_p(:));
        A_multi(:, pp) = utils.steering_vec_uca(M, radius_cal, lambda, ...
                            theta_cal_deg, phi_cal_multi_deg(pp));
    end
    [C_md_iR, c_md_iR, alpha_md, ~, n_iter_md] = ...
        estimate_C_multidir(B_multi, A_multi, M, 30, 1e-10);
    C_kw_multi{iRadius}  = C_md_iR;
    D_kw_multi{iRadius}  = inv(C_md_iR);
    err_C_kwmd(iRadius)  = norm(C_md_iR - C_true, 'fro') / norm(C_true, 'fro');
    c_md_vec{iRadius}    = c_md_iR;

    fprintf('  alpha_hat (KW 1-dir) = %+.4f %+.4fj\n', real(alpha_hat_iR), imag(alpha_hat_iR));
    fprintf('  ||residuo|| LS       = %.3e\n', norm(residual_iR));
    fprintf('  Erros Frobenius relativos:\n');
    fprintf('    Modelo IDEAL (oracle)        : %.4e   (sempre 0)\n', 0);
    fprintf('    Modelo perturb. %.0f%%          : %.4e\n', 100*pert_level_rel, err_Cb_F(iRadius));
    fprintf('    KW LS 1-dir                  : %.4e\n', err_C_kw1(iRadius));
    fprintf('    KW LS varias-dir P=%d (%d it.) : %.4e\n', P_multi, n_iter_md, err_C_kwmd(iRadius));
    fprintf('  (Coupling 6..9 sao SELF-CAL: estimam C dentro do loop principal)\n');

    fprintf('  Coeficientes c_k:\n');
    fprintf('    k |       true             |  KW LS 1-dir           |  KW LS varias-dir      |  modelo perturb.\n');
    fprintf('    --+------------------------+------------------------+------------------------+------------------------\n');
    for kk = 1:numel(c_true_vec)
        fprintf('    %d | %+.4f %+.4fj | %+.4f %+.4fj | %+.4f %+.4fj | %+.4f %+.4fj\n', kk, ...
            real(c_true_vec(kk)),  imag(c_true_vec(kk)), ...
            real(c_hat_iR(kk)),    imag(c_hat_iR(kk)), ...
            real(c_md_iR(kk)),     imag(c_md_iR(kk)), ...
            real(c_pert(kk)),      imag(c_pert(kk)));
    end
end

% --- Plot de diagnostico da calibracao ---
fig_cal = figure('Name', 'Calibracao - C_true vs estimativas (offline)', ...
                 'Color', 'w', 'Position', [50 50 1500 800]);

subplot(2,3,1);
imagesc(abs(C_true)); colorbar; axis equal tight;
title('|C_{true}|'); xlabel('coluna'); ylabel('linha');

subplot(2,3,2);
imagesc(abs(C_kw_1dir{1})); colorbar; axis equal tight;
title(sprintf('|C^{KW LS 1-dir}|  (err=%.2e)', err_C_kw1(1)));

subplot(2,3,3);
imagesc(abs(C_kw_multi{1})); colorbar; axis equal tight;
title(sprintf('|C^{KW LS Multi P=%d}|  (err=%.2e)', P_multi, err_C_kwmd(1)));

subplot(2,3,4);
imagesc(abs(C_model_b{1})); colorbar; axis equal tight;
title(sprintf('|C^{Modelo Pert. %.0f%%}|  (err=%.2e)', ...
    100*pert_level_rel, err_Cb_F(1)));
xlabel('coluna'); ylabel('linha');

subplot(2,3,5);
imagesc(abs(C_kw_1dir{1} - C_true)); colorbar; axis equal tight;
title('|C^{KW 1-dir} - C_{true}|');

subplot(2,3,6);
% Comparativo de erros (barras log)
err_data = [0, err_Cb_F(1), err_C_kw1(1), err_C_kwmd(1)];
labels_err = {'Ideal','Pert.','KW 1-dir', sprintf('KW Multi P=%d', P_multi)};
bar(err_data); set(gca, 'XTickLabel', labels_err, 'XTickLabelRotation', 15);
ylabel('||C - C_{true}||_F / ||C_{true}||_F');
set(gca,'YScale','log'); grid on;
title('Erros Frobenius relativos (offline)');

sgtitle(sprintf(['Calibracao OFFLINE (Coupling 2..5): SNR_{cal}=%d dB, ' ...
                 '\\phi_{cal}=%d°, N=%d  |  pert. modelo=%.0f%%  |  ' ...
                 'Multi: \\phi=[%s]°'], ...
    SNR_cal_dB, phi_cal_deg, N, 100*pert_level_rel, ...
    num2str(phi_cal_multi_deg)), 'FontWeight','bold');
exportgraphics(fig_cal, fullfile(outDir, 'calibracao_C_hat_vs_C_true.png'), 'Resolution', 200);

%% ======= DoA KW vs Delay and Sum vs Capon =======

phi_scan = -10:0.5:10;    % graus
BER_scan_in = zeros(size(phi_scan));
BER_scan_mvdr = zeros(size(phi_scan));
BER_scan_das = zeros(size(phi_scan));
EVM_scan_in = zeros(size(phi_scan));
EVM_scan_mvdr = zeros(size(phi_scan));
EVM_scan_das = zeros(size(phi_scan));

% Armazenamento dos resultados
phi_bp_all = cell(length(range_radius),1);
B_dB_DAS_all = cell(length(range_radius),1);
B_dB_Capon_all = cell(length(range_radius),1);
legend_entries = cell(length(range_radius),1);

% --- Estrutura de comparacao entre os cenarios ---
% Salva resultados no caso de referencia (phi_sig=75, phi_int=100)
results_cmp(nCoupling) = struct( ...
    'label',[], ...
    'phi_KW',NaN, 'phi_DAS',NaN, 'phi_MVDR',NaN, 'phi_MUSIC',NaN, ...
    'P_DAS_dB',[], 'P_MVDR_dB',[], 'P_MUSIC_dB',[], 'phi_scan',[], ...
    'BP_DAS_dB',[], 'BP_Capon_dB',[], 'phi_bp',[], ...
    'BER_in',NaN, 'BER_DAS',NaN, 'BER_MVDR',NaN, ...
    'EVM_in',NaN, 'EVM_DAS',NaN, 'EVM_MVDR',NaN, ...
    'sc_history',[]);   % NOVO: history struct retornado por estimate_C_selfcal
labels_cmp = ["No coupling (ideal)", ...
              "Coupling, no comp.", ...
              "Coupling, comp. modelo ideal", ...
              "Coupling, comp. modelo perturb.", ...
              "Coupling, comp. KW LS 1-dir", ...
              "Coupling, comp. KW LS varias-dir", ...
              "Coupling, self-cal DAS", ...
              "Coupling, self-cal CAPON", ...
              "Coupling, self-cal MUSIC", ...
              "Coupling, self-cal KW"];
for ic = 1:nCoupling
    results_cmp(ic).label = labels_cmp(range_coupling(ic)+1);
end

for iSNR = 1:nSNR
    for iISR = 1:nISR
        for iCoupling = 1:nCoupling
            for iRadius = 1:nRadius

                SNR_dB = range_SNR_dB(iSNR);
                ISR_dB = range_ISR_dB(iISR);
                radius = range_radius(iRadius)*lambda;
                Coupling = range_coupling(iCoupling);

                ii = 1;
                for iPhi = 1:nPhi
                    phi_sig_deg = range_phi(iPhi);

                    for iiPhi = 1:nPhi
                        phi_int_deg = range_phi(iiPhi);

                        if(phi_int_deg ~= phi_sig_deg)
                            j=1;

                            %phi_sig_deg = -180 + 360*rand;
                            %phi_int_deg = -180 + 360*rand;

                            % ----- Forma de onda conhecida (q) e interferidor "ruído" (r) -----

                            [X, q, r_int, Xsig, Xint, Xn, bits, pam_rrc_tx, pam_rect, sym_tx, qn, taus_sig, taus_int] = ...
                                utils.simulate_fsk_data_uca(M, radius, lambda, phi_sig_deg, phi_int_deg, ...
                                theta_sig_deg, theta_int_deg, SNR_dB, ISR_dB, N, fs, Rs, sps, alpha, span, fd);

                            % Aplicando o Mutual Coupling no canal
                            if Coupling == 0
                                % Caso 0: ideal, sem acoplamento
                                Coupling_matrix = eye(M);
                            else
                                % Casos 1..4: canal SEMPRE com acoplamento real (C_true).
                                % O que muda entre eles e' a estrategia de
                                % compensacao aplicada APOS a recepcao.
                                Coupling_matrix = Coupling_matrices(:,:,iRadius);
                            end

                            % X = Coupling_matrix * X;
                            % X = Coupling_matrix_sig * Xsig + Coupling_matrix_int * Xint + Xn;
                            X_1antenna = Xsig + Xint + Xn;
                            X_1antenna = X_1antenna(1,:);
                            Xq = Xsig + Xint + Xn;
                            % X = Coupling_matrix * Xsig + Coupling_matrix * Xint + Xn;
                            X = Coupling_matrix * Xsig + Coupling_matrix * Xint;
                            X = X ./ sqrt(mean(abs(X).^2));                           
                            X = X + Xn;

                            % --- Casos com compensacao (Coupling 2..9) ---
                            % O receptor compensa o acoplamento (sem saber a DoA do sinal).
                            sc_hist_iter = [];   % history p/ logging (so para self-cal)
                            switch Coupling
                                case 2   % modelo IDEAL (oracle)
                                    X = D_model_a{iRadius} * X;
                                case 3   % modelo PERTURBADO
                                    X = D_model_b{iRadius} * X;
                                case 4   % KW LS 1-dir
                                    X = D_kw_1dir{iRadius} * X;
                                case 5   % KW LS varias-dir
                                    X = D_kw_multi{iRadius} * X;
                                case {6, 7, 8, 9}   % SELF-CAL com DAS/CAPON/MUSIC/KW
                                    radius_m_sc = range_radius(iRadius)*lambda;
                                    sc_method = selfcal_methods{Coupling - 5};
                                    [C_sc, ~, ~, ~, sc_hist_iter] = estimate_C_selfcal(...
                                        X, q, M, radius_m_sc, lambda, ...
                                        sc_method, selfcal_grid_deg, ...
                                        selfcal_max_iter, [], [], C_true, ...
                                        selfcal_damping, selfcal_init);
                                    X = (C_sc \ X);   % equivalente a inv(C_sc)*X
                                % case 0 ou 1: nao faz nada
                            end

                            K = range_snapshots(end);
                            iTotal = iTotal + 1;
                            % ----- DoA KW
                            fprintf('--> Sim %d / %d , ii = %d , j = %d , K = %d / %d \n', iTotal, TotalSim, ...
                                ii, j, K, N_DOA);
                            % theta_hat_deg   = doa_kw_2007(X(:,1:K), q(1:K), M, d, lambda, 1);        % ell = 1 com M=4

                            beta = 2*pi*(0:M-1)'/M;
                            [theta_hat_deg, phi_hat_deg] = doa_kw_uca(X(:,1:K), q(1:K).', radius, lambda, beta);
                            % [theta_new_hat_deg, phi_new_hat_deg] = new_doa_kw_uca(X(:,1:K), q(1:K).', radius, lambda, beta);

                            fprintf('Estimativa KW: %+5.4f° | Vdd: %+5.4f°\n', phi_hat_deg, phi_sig_deg);

                            % Matriz de covariância
                            Rxx = (X(:,1:K)*X(:,1:K)')/size(X(:,1:K),2);
                            delta = 1e-3 * trace(Rxx)/M;
                            Rxx_dl = Rxx + delta*eye(M);

                            % Decomposição espectral
                            [eigvec, eigval] = eig(Rxx_dl);
                            [lambda_eig, idx] = sort(diag(eigval), 'descend');
                            E = eigvec(:, idx);

                            % Número de fontes conhecidas (K)
                            Ksrc = 2;   % ajuste se tiver mais fontes
                            En = E(:, Ksrc+1:end);  % subespaço do ruído

                            % Varredura angular (azimute, plano horizontal)
                            phi_scan = -180:0.1:180;    % graus
                            theta_scan = 90;            % fixa em 90° (plano XY)

                            P_DAS   = zeros(size(phi_scan));
                            P_MVDR  = zeros(size(phi_scan));
                            P_MUSIC = zeros(size(phi_scan));

                            Rinv = inv(Rxx_dl);

                            for ang = 1:numel(phi_scan)
                                % Estimando apenas o angulo phi
                                a = utils.steering_vec_uca(M, radius, lambda, theta_sig_deg, phi_scan(ang));  % Mx1
                                % ----- Delay-and-Sum -----
                                P_DAS(ang)  = abs(a' * Rxx * a);
                                % ----- Capon (MPDR) -----
                                denom     = real(a' * Rinv * a);
                                P_MVDR(ang) = 1 ./ max(denom, eps);
                                % ----- MUSIC -----
                                P_MUSIC(ang) = 1 ./ real(a' * (En * En') * a);
                            end

                            % Normalização (dB) e estimativa dos picos
                            P_DAS_dB   = 10*log10(P_DAS / max(P_DAS));
                            P_MVDR_dB  = 10*log10(P_MVDR / max(P_MVDR));
                            P_MUSIC_dB = 10*log10(P_MUSIC / max(P_MUSIC));

                            [~,i_das ] = max(P_DAS);
                            [~,i_mvdr] = max(P_MVDR);
                            [~, idx_music] = max(P_MUSIC);

                            phi_DAS  = phi_scan(i_das);
                            phi_MVDR = phi_scan(i_mvdr);
                            phi_MUSIC = phi_scan(idx_music);

                            % --- NOVO: armazenar para comparacao final
                            % (caso de referencia: phi_sig=75) ---
                            if abs(phi_sig_deg - 75) < 1e-9
                                results_cmp(iCoupling).phi_KW    = phi_hat_deg;
                                results_cmp(iCoupling).phi_DAS   = phi_DAS;
                                results_cmp(iCoupling).phi_MVDR  = phi_MVDR;
                                results_cmp(iCoupling).phi_MUSIC = phi_MUSIC;
                                results_cmp(iCoupling).P_DAS_dB   = P_DAS_dB;
                                results_cmp(iCoupling).P_MVDR_dB  = P_MVDR_dB;
                                results_cmp(iCoupling).P_MUSIC_dB = P_MUSIC_dB;
                                results_cmp(iCoupling).phi_scan   = phi_scan;
                                % Para self-cal (Coupling 6..9), salva history
                                if ~isempty(sc_hist_iter)
                                    results_cmp(iCoupling).sc_history = sc_hist_iter;
                                end
                            end

                            switch Coupling
                                case 0, name_string = 'No Coupling';
                                case 1, name_string = 'Coupling (uncomp.)';
                                case 2, name_string = 'Coupling (comp. ideal)';
                                case 3, name_string = 'Coupling (comp. pert.)';
                                case 4, name_string = 'Coupling (KW LS 1-dir)';
                                case 5, name_string = 'Coupling (KW LS Multi)';
                                case 6, name_string = 'Coupling (Self-cal DAS)';
                                case 7, name_string = 'Coupling (Self-cal CAPON)';
                                case 8, name_string = 'Coupling (Self-cal MUSIC)';
                                case 9, name_string = 'Coupling (Self-cal KW)';
                            end
                             % ----- PLOT DOA -----
                            % fig_name = "DoA DAS (" + string(phi_sig_deg) + " graus) with " + name_string + ...
                            %     " (SNR=" + string(range_SNR_dB(iSNR)) + " | ISR=" + string(range_ISR_dB(iISR)) + " | r=" + string(range_radius(iRadius)) + ")";
                            % fig = figure( 'Name', fig_name, 'NumberTitle', 'off');
                            % hold on; grid on;
                            % 
                            % plot(phi_scan, P_MVDR_dB, 'LineWidth', 2);
                            % plot(phi_scan, P_DAS_dB, 'LineWidth', 2);
                            % xline(phi_sig_deg, '--r', sprintf('\\phi=%.1f°', phi_sig_deg));
                            % legend('Capon/MVDR', 'DAS','Location','best');
                            % % fileName = fig_name + ".png";
                            % % % Limpar nome
                            % % fileName = replace(fileName, "|", "-");
                            % % fileName = regexprep(fileName, '[^a-zA-Z0-9 _\-\.\(\)]', '');
                            % % exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);

                            % ----- ESTIMACAO DA MCM (DESATIVADO) -----
                            % O bloco abaixo (com fminsearch) foi substituido pela
                            % CALIBRACAO LS LINEAR feita uma unica vez no inicio do
                            % script (ver bloco "CALIBRACAO" antes do loop).
                            %
                            % O artigo de Khan 2020 mostra que UMA medicao em uma
                            % direcao conhecida ja basta para estimar D, pois o
                            % acoplamento depende SO da geometria. Logo, nao faz
                            % sentido re-estimar dentro do loop.
                            %
                            % if abs(phi_sig_deg - 75) < 1e-9 && Coupling == 1
                            %     a_sig = utils.steering_vec_uca(M, radius, lambda, theta_sig_deg, phi_sig_deg);
                            %     b_hat = X * q / (q' * q);
                            %     p0 = zeros(1,8);
                            %     cost_fun = @(p) cost_mcm_circulant4(p, a, b_hat);
                            %     opts = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10,...
                            %         'MaxIter',5000,'MaxFunEvals',20000);
                            %     p_est = fminsearch(cost_fun, p0, opts);
                            %     c1_est = p_est(1) + 1j*p_est(2);
                            %     c2_est = p_est(3) + 1j*p_est(4);
                            %     c3_est = p_est(5) + 1j*p_est(6);
                            %     c4_est = p_est(7) + 1j*p_est(8);
                            %     C_est = build_C_circulant_uca4(c1_est, c2_est, c3_est, c4_est);
                            %     fprintf('Erro relativo Frobenius (in-loop) = %.6e\n', ...
                            %         norm(C_est - C_true, 'fro') / norm(C_true, 'fro'));
                            % end

                            % ----- BEAMFORMING -----

                            % Varredura angular (azimute, plano horizontal)
                            phi_scan = -180:0.5:180;    % graus
                            % phi_scan = phi_hat_deg-30:0.5:phi_hat_deg+30;
                            for ang = 1:numel(phi_scan)
                                % ----- Restrições angulares
                                a_sig = utils.steering_vec_uca(M, radius, lambda, theta_sig_deg, phi_scan(ang));  % vetor de direção Mx1
                                % --- Direções dos interferentes (onde queremos forçar nulos) ---
                                phi_nulls = [phi_int_deg];                 % vetor com 1 ou mais ângulos de nulo (ex: [100 120])
                                Nnull = length(phi_nulls);

                                % --- Matriz de restrições C ---
                                % Cada coluna de C é um vetor de direção (sinal + interferentes)
                                C = a_sig;                                   % primeira coluna = sinal desejado
                                if nulo_no_interferidor
                                    for k = 1:Nnull
                                        a_null = utils.steering_vec_uca(M, radius, lambda, theta_sig_deg, phi_nulls(k));  % vetor de direção Mx1
                                        C = [C, a_null];                          % adiciona interferente
                                    end
                                end

                                % --- Vetor de ganhos desejados f ---
                                % Ganho = 1 para o sinal desejado, 0 para cada nulo
                                if nulo_no_interferidor
                                    f = [1; zeros(Nnull,1)];
                                else
                                    f = [1];
                                end

                                % --- Matriz de projeção e vetor de correção
                                P   = eye(M) - C * ((C' * C) \ C');
                                f_c = C * ((C' * C) \ f);

                                % ---- Capon ----
                                Rxx = (X(:,1:N_DOA)*X(:,1:N_DOA)')/N_DOA;
                                delta = 1e-3 * trace(Rxx)/M;
                                Rinv = inv(Rxx + delta*eye(M));
                                w_capon = (Rinv*a_sig) / (a_sig' * Rinv * a_sig);
                                y_capon = w_capon' * X;

                                % ---- Delay-and-Sum (DAS) ----
                                w_das = a_sig / M;
                                y_das = w_das' * X;       % saída DAS

                                % ======= Beampattern em 90° para conhecimento =======
                                if abs(phi_scan(ang) - teste_phi) < 1e-9 && abs(phi_sig_deg - 75) < 1e-9
                                    [phi_beampattern, B_dB_DAS] = utils.beampattern_db_uca(w_das, M, radius, lambda, phi_scan);
                                    [~,               B_dB_Capon] = utils.beampattern_db_uca(w_capon, M, radius, lambda, phi_scan);

                                    phi_bp_all{iRadius}      = phi_beampattern;
                                    B_dB_DAS_all{iRadius}    = B_dB_DAS;
                                    B_dB_Capon_all{iRadius}  = B_dB_Capon;

                                    legend_entries{iRadius} = sprintf('raio = %.2f m', range_radius(iRadius));

                                    % --- NOVO: salvar para comparacao final ---
                                    results_cmp(iCoupling).phi_bp      = phi_beampattern;
                                    results_cmp(iCoupling).BP_DAS_dB   = B_dB_DAS;
                                    results_cmp(iCoupling).BP_Capon_dB = B_dB_Capon;
                                end

                                % ======= PLOT DIAGNOSTICO: Sinal Antes/Depois do Acoplamento + Beamforming =======
                                if abs(phi_scan(ang) - teste_phi) < 1e-9 && abs(phi_sig_deg - 75) < 1e-9
                                    % --- Configuracao do plot ---
                                    ant_idx    = 1;           % indice da antena a plotar
                                    N_plot     = 1500;        % numero de amostras a mostrar (trecho inicial)

                                    % --- Sinais ---
                                    sig_ideal    = real(Xq(ant_idx, 1:N_plot));      % sem acoplamento (antena ant_idx)
                                    sig_acoplado = real(X(ant_idx, 1:N_plot));        % com acoplamento (antena ant_idx)
                                    sig_capon    = real(y_capon(1:N_plot));            % saida Capon (escalar)
                                    sig_das      = real(y_das(1:N_plot));              % saida DAS (escalar)

                                    % ---- Figura: IQ ----
                                    if Coupling == 1
                                    fig_name = 'Cadeia de Sinal: Ideal -> Acoplado -> Beamforming'
                                    else 
                                    fig_name = 'Cadeia de Sinal: Ideal -> Beamforming'
                                    end
                                    % figure('Position', [50 50 1500 900], 'Color', 'w', ...
                                    %     'Name', fig_name);
                                    % 
                                    % % (1) Sinal ideal (sem acoplamento)
                                    % subplot(5,1,1);
                                    % plot(real(q(1:N_plot)), 'b-', 'LineWidth', 0.8); hold on;
                                    % plot(imag(q(1:N_plot)), 'r-', 'LineWidth', 0.8);
                                    % grid on;
                                    % ylabel('Amplitude');
                                    % title(sprintf('(a) Sem ideal', ant_idx));
                                    % % legend('Re\{x_{ideal}(t)\}', 'Simbolo', 'Location', 'northeast');
                                    % 
                                    % % (1) Sinal ideal (sem acoplamento)
                                    % subplot(5,1,2);
                                    % plot(real(Xq(ant_idx,1:N_plot)), 'b-', 'LineWidth', 0.8); hold on;
                                    % plot(imag(Xq(ant_idx,1:N_plot)), 'r-', 'LineWidth', 0.8);
                                    % grid on;
                                    % ylabel('Amplitude');
                                    % title(sprintf('(b) Sem acoplamento - Antena %d', ant_idx));
                                    % % legend('Re\{x_{ideal}(t)\}', 'Simbolo', 'Location', 'northeast');
                                    % 
                                    % % (2) Sinal com acoplamento (sem beamforming)
                                    % subplot(5,1,3);
                                    % plot(real(X(ant_idx,1:N_plot)), 'b-', 'LineWidth', 0.8); hold on;
                                    % plot(imag(X(ant_idx,1:N_plot)), 'r-', 'LineWidth', 0.8);
                                    % grid on;
                                    % ylabel('Amplitude');
                                    % title(sprintf('(b) Com acoplamento (C^{-1}) - Antena %d', ant_idx));
                                    % % legend('Re\{x_{acoplado}(t)\}', 'Simbolo', 'Location', 'northeast');
                                    % 
                                    % % (3) Saida Capon (MVDR)
                                    % subplot(5,1,4);
                                    % plot(real(y_capon(1:N_plot)), 'b-', 'LineWidth', 0.8); hold on;
                                    % plot(imag(y_capon(1:N_plot)), 'r-', 'LineWidth', 0.8);
                                    % grid on;
                                    % ylabel('Amplitude');
                                    % title('(c) Saida Capon (MVDR) - apos acoplamento + beamforming');
                                    % % legend('Re\{y_{Capon}(t)\}', 'Simbolo', 'Location', 'northeast');
                                    % 
                                    % % (4) Saida DAS
                                    % subplot(5,1,5);
                                    % plot(real(y_das(1:N_plot)), 'b-', 'LineWidth', 0.8); hold on;
                                    % plot(imag(y_das(1:N_plot)), 'r-', 'LineWidth', 0.8);
                                    % grid on;
                                    % hold off; grid on;
                                    % xlabel('Tempo (ms)');
                                    % ylabel('Amplitude');
                                    % title('(d) Saida DAS - apos acoplamento + beamforming');
                                    % % legend('Re\{y_{DAS}(t)\}', 'Simbolo', 'Location', 'northeast');
                                    % 
                                    % sgtitle(sprintf('Cadeia de sinal | Coupling %d | SNR=%d dB | ISR=%d dB | R=%.2f\\lambda | \\phi_{sig}=%.0f° | \\phi_{int}=%.0f°', ...
                                    %     Coupling, SNR_dB, ISR_dB, range_radius(iRadius), phi_sig_deg, phi_int_deg), ...
                                    %     'FontSize', 13, 'FontWeight', 'bold');

                                    % ---- Figura 2: Sobreposicao de todos ----

                                end

                                % ======= Demodulação =======
                                [bits_hat_in, BER_in, pam_rx_mf1, sym_rx1]    = utils.fsk2_demod(X(1,:), bits, Rs, sps, alpha, span, fd);
                                [bits_hat_das, BER_das, pam_rx_mf2, sym_rx2]  = utils.fsk2_demod(y_das, bits, Rs, sps, alpha, span, fd);
                                [bits_hat_mvdr, BER_mvdr, pam_rx_mf3, sym_rx3]= utils.fsk2_demod(y_capon, bits, Rs, sps, alpha, span, fd);

                                % Calcula BER
                                BER_scan_in(ang) = BER_in*100;
                                BER_scan_mvdr(ang) = BER_mvdr*100;
                                BER_scan_das(ang) = BER_das*100;

                                fprintf('BER: RX %.2f%% / Capon %.2f%% / DAS %.2f%% \n', ...
                                    BER_in*100, BER_mvdr*100, BER_das*100);

                                % Calcula EVM
                                [EVM_in,  EVMdB_in]                 = utils.calc_evm_real(sym_rx1,  sym_tx);
                                [EVM_das, EVMdB_das]                = utils.calc_evm_real(sym_rx2, sym_tx);
                                [EVM_mvdr,EVMdB_mvdr]               = utils.calc_evm_real(sym_rx3, sym_tx);

                                EVM_scan_in(ang) = 20*log10(EVM_in);
                                EVM_scan_mvdr(ang) = 20*log10(EVM_mvdr);
                                EVM_scan_das(ang) = 20*log10(EVM_das);

                                fprintf('EVM (dB): Rx: %.2f | Capon: %.2f | DAS: %.2f\n', ...
                                    EVMdB_in, EVMdB_mvdr, EVMdB_das);

                                % --- NOVO: salvar BER/EVM para comparacao final entre cenarios ---
                                if abs(phi_scan(ang) - teste_phi) < 1e-9 && abs(phi_sig_deg - 75) < 1e-9
                                    results_cmp(iCoupling).BER_in   = BER_in   * 100;
                                    results_cmp(iCoupling).BER_DAS  = BER_das  * 100;
                                    results_cmp(iCoupling).BER_MVDR = BER_mvdr * 100;
                                    results_cmp(iCoupling).EVM_in   = 20*log10(EVM_in);
                                    results_cmp(iCoupling).EVM_DAS  = 20*log10(EVM_das);
                                    results_cmp(iCoupling).EVM_MVDR = 20*log10(EVM_mvdr);
                                end

                            end


                            if abs(phi_sig_deg - 75) < 1e-9
                                fig_name = "BER, EVM (" + string(phi_sig_deg) + "°) " + name_string + ...
                                    " (r=" + string(range_radius(iRadius))  + " | ISR=" + string(range_ISR_dB(iISR)) + " | SNR=" + string(range_SNR_dB(iSNR)) + ")";
                                fig = figure( 'Name', fig_name, 'NumberTitle', 'off');

                                % ======= SUBPLOT 1: BER =======
                                subplot(4,1,1);  % 2 linhas, 1 coluna, posição 1
                                hold on; grid on;

                                plot(phi_scan, BER_scan_in, 'Color', colors(1), 'LineWidth', 2);
                                plot(phi_scan, BER_scan_mvdr, 'Color', colors(2), 'LineWidth', 2);
                                plot(phi_scan, BER_scan_das, 'Color', colors(3), 'LineWidth', 2);

                                [~, idx_mvdr] = min(BER_scan_mvdr);
                                [~, idx_das]  = min(BER_scan_das);

                                xline(phi_MVDR, '--g', sprintf('\\phi_{MVDR}=%.1f°', phi_MVDR));
                                xline(phi_DAS,  '--b', sprintf('\\phi_{DAS}=%.1f°', phi_DAS));
                                xline(phi_hat_deg, '--k', sprintf('\\phi_{KW}=%.1f°', phi_hat_deg));
                                xline(phi_sig_deg, '--r', sprintf('\\phi_{real}=%.1f°', phi_sig_deg));

                                ylabel("BER (%)",'Interpreter','latex');
                                legend(methodsBeamforming, 'Location', 'northwest','Interpreter','latex');
                                title('BER');

                                ylim([0 55]);
                                % xlim([-180 180]);

                                % ======= SUBPLOT 2: EVM =======
                                subplot(4,1,2);  % posição 2
                                hold on; grid on;

                                plot(phi_scan, EVM_scan_in, 'Color', colors(1), 'LineWidth', 2);
                                plot(phi_scan, EVM_scan_mvdr, 'Color', colors(2), 'LineWidth', 2);
                                plot(phi_scan, EVM_scan_das, 'Color', colors(3), 'LineWidth', 2);

                                [~, idx_mvdr] = min(EVM_scan_mvdr);
                                [~, idx_das]  = min(EVM_scan_das);

                                xline(phi_MVDR, '--g', sprintf('\\phi_{MVDR}=%.1f°', phi_MVDR));
                                xline(phi_DAS,  '--b', sprintf('\\phi_{DAS}=%.1f°', phi_DAS));
                                xline(phi_hat_deg, '--k', sprintf('\\phi_{KW}=%.1f°', phi_hat_deg));
                                xline(phi_sig_deg, '--r', sprintf('\\phi_{real}=%.1f°', phi_sig_deg));

                                xlabel("angulo (graus)",'Interpreter','latex');
                                ylabel("EVM (dB)",'Interpreter','latex');
                                legend(methodsBeamforming, 'Location', 'northwest','Interpreter','latex');
                                title('EVM');

                                switch Coupling
                                    case 0, name_string = 'No Coupling';
                                    case 1, name_string = 'Coupling (uncomp.)';
                                    case 2, name_string = 'Coupling (comp. ideal)';
                                    case 3, name_string = 'Coupling (comp. pert.)';
                                    case 4, name_string = 'Coupling (KW LS 1-dir)';
                                    case 5, name_string = 'Coupling (KW LS Multi)';
                                    case 6, name_string = 'Coupling (Self-cal DAS)';
                                    case 7, name_string = 'Coupling (Self-cal CAPON)';
                                    case 8, name_string = 'Coupling (Self-cal MUSIC)';
                                    case 9, name_string = 'Coupling (Self-cal KW)';
                                end

                                % ======= SUBPLOT 3: BEAMPATTERN =======
                                subplot(4,1,3);  % posição 2
                                hold on; grid on;

                                for ir = 1:length(range_radius)
                                    if ~isempty(B_dB_DAS_all{ir})
                                        plot(phi_bp_all{ir}, B_dB_DAS_all{ir}, 'LineWidth', 1.4);
                                    end
                                end

                                xline(75, 'g--', 'phi\_hat\_deg');
                                xline(teste_phi, 'g--', 'Teste');
                                xlabel('\phi (graus)');
                                ylabel('Resposta (dB)');
                                title('Padrão de radiação - DAS');
                                legend(legend_entries(~cellfun('isempty', legend_entries)), 'Location', 'best');

                                ylim([-40 5]);
                                xlim([-180 180]);

                                subplot(4,1,4);  % posição 2
                                hold on; grid on;
                                for ir = 1:length(range_radius)
                                    if ~isempty(B_dB_Capon_all{ir})
                                        plot(phi_bp_all{ir}, B_dB_Capon_all{ir}, 'LineWidth', 1.4);
                                    end
                                end

                                xline(75, 'g--', 'phi\_hat\_deg');
                                xline(teste_phi, 'g--', 'Teste');
                                xlabel('\phi (graus)');
                                ylabel('Resposta (dB)');
                                title('Padrão de radiação - Capon');
                                legend(legend_entries(~cellfun('isempty', legend_entries)), 'Location', 'best');

                                ylim([-40 5]);
                                xlim([-180 180]);


                                fileName = fig_name + ".png";
                                % Limpar nome
                                fileName = replace(fileName, "|", "-");
                                fileName = regexprep(fileName, '[^a-zA-Z0-9 _\-\.\(\)]', '');

                                exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);
                            end
                            % Salvando EVM e BER
                            % BER(:,iSNR,iISR,iRadius, iCoupling, m, ii) = [
                            %     BER_in;
                            %     BER_das
                            %     ];
                            %
                            % EVM(:,iSNR,iISR,iRadius, iCoupling, m, ii) = [
                            %     EVMdB_in;
                            %     EVMdB_das
                            %     ];
                        end
                        ii = ii + 1;

                    end
                end
            end
            % if Coupling == 1
            %     name_string = 'Coupling';
            % else
            %     name_string = 'No Coupling';
            % end
            % 
            % % ======= Plot BEAMPATTERN =======
            % 
            % figure('Name',"Beampattern DAS with " + name_string);
            % hold on; grid on; box on;
            % 
            % for ir = 1:length(range_radius)
            %     if ~isempty(B_dB_DAS_all{ir})
            %         plot(phi_bp_all{ir}, B_dB_DAS_all{ir}, 'LineWidth', 1.4);
            %     end
            % end
            % 
            % xline(75, 'g--', 'phi\_hat\_deg');
            % 
            % xlabel('\phi (graus)');
            % ylabel('Resposta (dB)');
            % title('Padrão de radiação - DAS');
            % legend(legend_entries(~cellfun('isempty', legend_entries)), 'Location', 'best');
            % 
            % ylim([-40 5]);
            % xlim([-180 180]);
            % 
            % % ---
            % figure('Name',"Beampattern Capon with " + name_string);
            % hold on; grid on; box on;
            % 
            % for ir = 1:length(range_radius)
            %     if ~isempty(B_dB_Capon_all{ir})
            %         plot(phi_bp_all{ir}, B_dB_Capon_all{ir}, 'LineWidth', 1.4);
            %     end
            % end
            % 
            % xline(75, 'g--', 'phi\_hat\_deg');
            % 
            % xlabel('\phi (graus)');
            % ylabel('Resposta (dB)');
            % title('Padrão de radiação - Capon');
            % legend(legend_entries(~cellfun('isempty', legend_entries)), 'Location', 'best');
            % 
            % ylim([-40 5]);
            % xlim([-180 180]);
        end
    end
end

% =========================================================================
% PLOT FINAL: Comparacao dos cenarios (sem acoplamento / com / compensado)
% phi_sig = 75 deg, phi_int = 100 deg, teste_phi = 75
% =========================================================================
phi_sig_ref = 75;
phi_int_ref = 100;
teste_phi_ref = teste_phi;

% --- Identifica quais cenarios foram efetivamente preenchidos ---
valid = false(nCoupling,1);
for ic = 1:nCoupling
    valid(ic) = ~isempty(results_cmp(ic).phi_scan) && ~isnan(results_cmp(ic).phi_KW);
end
ic_valid = find(valid);

if ~isempty(ic_valid)
    % --- Tabela de estimativas de DoA com erro absoluto ---
    fprintf('\n=========================================================\n');
    fprintf('Comparacao final entre cenarios (phi_sig=%g, phi_int=%g)\n', ...
            phi_sig_ref, phi_int_ref);
    fprintf('=========================================================\n');
    fprintf('%-25s | %8s | %8s | %8s | %8s\n', ...
            'Cenario','KW','DAS','MVDR','MUSIC');
    fprintf('---------------------------------------------------------\n');
    for k = 1:numel(ic_valid)
        ic = ic_valid(k);
        fprintf('%-25s | %8.3f | %8.3f | %8.3f | %8.3f\n', ...
            char(results_cmp(ic).label), ...
            results_cmp(ic).phi_KW, results_cmp(ic).phi_DAS, ...
            results_cmp(ic).phi_MVDR, results_cmp(ic).phi_MUSIC);
    end
    fprintf('---------------------------------------------------------\n');
    fprintf('Erro absoluto (graus) vs phi_sig = %g:\n', phi_sig_ref);
    fprintf('%-25s | %8s | %8s | %8s | %8s\n', ...
            'Cenario','|eKW|','|eDAS|','|eMVDR|','|eMUSIC|');
    fprintf('---------------------------------------------------------\n');
    for k = 1:numel(ic_valid)
        ic = ic_valid(k);
        fprintf('%-25s | %8.3f | %8.3f | %8.3f | %8.3f\n', ...
            char(results_cmp(ic).label), ...
            abs(results_cmp(ic).phi_KW    - phi_sig_ref), ...
            abs(results_cmp(ic).phi_DAS   - phi_sig_ref), ...
            abs(results_cmp(ic).phi_MVDR  - phi_sig_ref), ...
            abs(results_cmp(ic).phi_MUSIC - phi_sig_ref));
    end
    fprintf('=========================================================\n\n');

    % --- Figura comparativa ---
    fig_cmp = figure('Name','Comparacao final: cenarios de acoplamento', ...
                     'NumberTitle','off', 'Position',[50 50 1700 950]);

    cmap = lines(numel(ic_valid));

    % Labels curtos para os eixos das barras
    short_labels_map = containers.Map( ...
        {'No coupling (ideal)', ...
         'Coupling, no comp.', ...
         'Coupling, comp. modelo ideal', ...
         'Coupling, comp. modelo perturb.', ...
         'Coupling, comp. KW LS 1-dir', ...
         'Coupling, comp. KW LS varias-dir', ...
         'Coupling, self-cal DAS', ...
         'Coupling, self-cal CAPON', ...
         'Coupling, self-cal MUSIC', ...
         'Coupling, self-cal KW'}, ...
        {'No Coup.', 'No comp.', 'Ideal', 'Pert.', 'KW 1-dir', 'KW Multi', ...
         'SC-DAS', 'SC-CAPON', 'SC-MUSIC', 'SC-KW'});
    short_labels = cell(numel(ic_valid),1);
    for k = 1:numel(ic_valid)
        ic = ic_valid(k);
        key = char(results_cmp(ic).label);
        if isKey(short_labels_map, key)
            short_labels{k} = short_labels_map(key);
        else
            short_labels{k} = key;
        end
    end

    % (1) Espectros MUSIC sobrepostos
    subplot(2,2,1); hold on; grid on;
    for k = 1:numel(ic_valid)
        ic = ic_valid(k);
        plot(results_cmp(ic).phi_scan, results_cmp(ic).P_MUSIC_dB, ...
             'LineWidth', 1.6, 'Color', cmap(k,:), ...
             'DisplayName', char(results_cmp(ic).label));
    end
    xline(phi_sig_ref, 'k--', 'phi_{sig}', 'LabelVerticalAlignment','bottom');
    xline(phi_int_ref, 'r--', 'phi_{int}', 'LabelVerticalAlignment','bottom');
    xlabel('\phi (graus)'); ylabel('Pseudoespectro MUSIC (dB)');
    title('MUSIC: comparacao entre cenarios');
    legend('Location','best'); xlim([min(results_cmp(ic_valid(1)).phi_scan) ...
                                     max(results_cmp(ic_valid(1)).phi_scan)]);

    % (2) Beampatterns Capon sobrepostos
    subplot(2,2,2); hold on; grid on;
    for k = 1:numel(ic_valid)
        ic = ic_valid(k);
        if ~isempty(results_cmp(ic).BP_Capon_dB)
            plot(results_cmp(ic).phi_bp, results_cmp(ic).BP_Capon_dB, ...
                 'LineWidth', 1.6, 'Color', cmap(k,:), ...
                 'DisplayName', char(results_cmp(ic).label));
        end
    end
    xline(phi_sig_ref, 'k--', 'phi_{sig}', 'LabelVerticalAlignment','bottom');
    xline(phi_int_ref, 'r--', 'phi_{int}', 'LabelVerticalAlignment','bottom');
    xlabel('\phi (graus)'); ylabel('|w^H a(\phi)|^2 (dB)');
    title(sprintf('Beampattern Capon (apontado p/ %g°)', teste_phi_ref));
    legend('Location','best'); ylim([-60 5]);

    % (3) Barras BER
    subplot(2,2,3);
    BER_mat = zeros(numel(ic_valid), 3);
    for k = 1:numel(ic_valid)
        ic = ic_valid(k);
        BER_mat(k,:) = [results_cmp(ic).BER_in, ...
                        results_cmp(ic).BER_DAS, ...
                        results_cmp(ic).BER_MVDR];
    end
    bh = bar(BER_mat, 'grouped');
    set(gca, 'XTickLabel', short_labels, ...
        'XTickLabelRotation', 25);
    grid on;
    ylabel('BER (%)');
    title('BER por cenario × beamformer');
    legend({'Rx (1 antena)','DAS','MVDR/Capon'}, 'Location','best');

    % (4) Barras EVM
    subplot(2,2,4);
    EVM_mat = zeros(numel(ic_valid), 3);
    for k = 1:numel(ic_valid)
        ic = ic_valid(k);
        EVM_mat(k,:) = [results_cmp(ic).EVM_in, ...
                        results_cmp(ic).EVM_DAS, ...
                        results_cmp(ic).EVM_MVDR];
    end
    bh2 = bar(EVM_mat, 'grouped');
    set(gca, 'XTickLabel', short_labels, ...
        'XTickLabelRotation', 25);
    grid on;
    ylabel('EVM (dB)');
    title('EVM por cenario × beamformer');
    legend({'Rx (1 antena)','DAS','MVDR/Capon'}, 'Location','best');

    sgtitle(sprintf(['Comparacao: 10 cenarios de acoplamento\n' ...
                     '\\phi_{sig}=%g°, \\phi_{int}=%g°, SNR=%g dB, ISR=%g dB, ' ...
                     'raio=%.2f \\lambda, pert. modelo=%.0f%%, KW Multi P=%d'], ...
            phi_sig_ref, phi_int_ref, range_SNR_dB(end), range_ISR_dB(end), ...
            range_radius(end), 100*pert_level_rel, P_multi));

    exportgraphics(fig_cmp, fullfile(outDir, 'comparacao_cenarios_acoplamento.png'), ...
                   'Resolution', 200);

    % =====================================================================
    % PLOT EXTRA: Curva de aprendizado da self-cal (4 metodos de DoA)
    % Mostra erro Frobenius || C_t - C_true ||_F / ||C_true||_F vs iteracao,
    % alem da DoA estimada por iteracao e o salto |C_t - C_{t-1}|.
    % =====================================================================
    sc_indices_in_results = [];   % indices dentro de results_cmp que tem sc_history
    sc_method_names = {};
    for ic = 1:nCoupling
        coup_code = range_coupling(ic);
        if coup_code >= 6 && coup_code <= 9 && ~isempty(results_cmp(ic).sc_history)
            sc_indices_in_results(end+1) = ic; %#ok<AGROW>
            sc_method_names{end+1}       = selfcal_methods{coup_code - 5}; %#ok<AGROW>
        end
    end

    if ~isempty(sc_indices_in_results)
        fig_lc = figure('Name','Learning curve - Self-cal por metodo de DoA', ...
                        'NumberTitle','off', 'Position',[100 100 1400 850]);
        cmap_sc = lines(numel(sc_indices_in_results));

        % --- (1) Erro Frobenius vs iteracao -------------------------------
        subplot(2,2,1); hold on; grid on;
        for kk = 1:numel(sc_indices_in_results)
            ic = sc_indices_in_results(kk);
            h = results_cmp(ic).sc_history;
            it_axis = 1:numel(h.err_F_per_iter);
            plot(it_axis, h.err_F_per_iter, 'o-', ...
                 'LineWidth', 1.8, 'Color', cmap_sc(kk,:), ...
                 'MarkerFaceColor', cmap_sc(kk,:), ...
                 'DisplayName', sprintf('SC-%s', sc_method_names{kk}));
        end
        % Linha de referencia: erro do KW LS 1-dir
        kw1_idx = find(range_coupling == 4, 1);
        if ~isempty(kw1_idx)
            yline(err_C_kw1(1), 'k--', 'LineWidth', 1.2, ...
                  'DisplayName', 'KW LS 1-dir (offline)');
        end
        kwmd_idx = find(range_coupling == 5, 1);
        if ~isempty(kwmd_idx)
            yline(err_C_kwmd(1), 'k:', 'LineWidth', 1.2, ...
                  'DisplayName', sprintf('KW LS Multi P=%d (offline)', P_multi));
        end
        set(gca, 'YScale', 'log');
        xlabel('iteracao'); ylabel('||C_t - C_{true}||_F / ||C_{true}||_F');
        title('Curva de aprendizado: erro Frobenius vs iteracao');
        legend('Location','northeast');

        % --- (2) DoA estimada por iteracao --------------------------------
        subplot(2,2,2); hold on; grid on;
        for kk = 1:numel(sc_indices_in_results)
            ic = sc_indices_in_results(kk);
            h = results_cmp(ic).sc_history;
            it_axis = 1:numel(h.phi_per_iter);
            plot(it_axis, h.phi_per_iter, 'o-', ...
                 'LineWidth', 1.6, 'Color', cmap_sc(kk,:), ...
                 'MarkerFaceColor', cmap_sc(kk,:), ...
                 'DisplayName', sprintf('SC-%s', sc_method_names{kk}));
        end
        yline(phi_sig_ref, 'k--', 'LineWidth', 1.4, 'DisplayName','\phi_{sig}');
        xlabel('iteracao'); ylabel('\phi estimada (°)');
        title('DoA estimada por iteracao');
        legend('Location','best');

        % --- (3) |C_t - C_{t-1}| (criterio interno de convergencia) -------
        subplot(2,2,3); hold on; grid on;
        for kk = 1:numel(sc_indices_in_results)
            ic = sc_indices_in_results(kk);
            h = results_cmp(ic).sc_history;
            it_axis = 1:numel(h.delta_C);
            plot(it_axis, h.delta_C, 'o-', ...
                 'LineWidth', 1.6, 'Color', cmap_sc(kk,:), ...
                 'MarkerFaceColor', cmap_sc(kk,:), ...
                 'DisplayName', sprintf('SC-%s', sc_method_names{kk}));
        end
        set(gca, 'YScale', 'log');
        xlabel('iteracao'); ylabel('||C_t - C_{t-1}||_F / ||C_{t-1}||_F');
        title('Salto entre iteracoes (criterio de parada)');
        legend('Location','best');

        % --- (4) Tabela de iteracoes ate convergencia ---------------------
        subplot(2,2,4); axis off;
        txt = {'\bf{Resumo Self-Cal:}', ''};
        txt{end+1} = sprintf('  %-10s | %5s | %s', 'metodo', 'n_it', 'err.Frob.final');
        txt{end+1} = repmat('-', 1, 45);
        for kk = 1:numel(sc_indices_in_results)
            ic = sc_indices_in_results(kk);
            h  = results_cmp(ic).sc_history;
            txt{end+1} = sprintf('  %-10s | %5d | %.3e', ...
                sprintf('SC-%s', sc_method_names{kk}), h.n_iter, ...
                h.err_F_per_iter(end));
        end
        txt{end+1} = '';
        txt{end+1} = '\bf{Referencias offline:}';
        txt{end+1} = sprintf('  KW LS 1-dir : %.3e', err_C_kw1(1));
        txt{end+1} = sprintf('  KW LS Multi : %.3e', err_C_kwmd(1));
        txt{end+1} = sprintf('  Modelo pert.: %.3e', err_Cb_F(1));
        text(0.05, 0.95, txt, 'Units','normalized', ...
             'VerticalAlignment','top', 'FontName','Courier', ...
             'Interpreter','tex', 'FontSize', 11);

        sgtitle(sprintf(['Self-cal: convergencia por metodo de DoA\n' ...
                         '\\phi_{sig}=%g°, \\phi_{int}=%g°, SNR=%g dB, ISR=%g dB'], ...
                phi_sig_ref, phi_int_ref, range_SNR_dB(end), range_ISR_dB(end)), ...
                'FontWeight','bold');

        exportgraphics(fig_lc, fullfile(outDir, 'selfcal_learning_curve.png'), ...
                       'Resolution', 200);
    end
else
    warning(['Nenhum cenario foi preenchido em results_cmp. ' ...
             'Verifique se phi_sig=75 e teste_phi=75 estao no range simulado.']);
end

return;

RMSE_mean = mean(RMSE, 7);   % média ao longo da 7ª dimensão (trials)
BER_mean = mean(BER, 7);   % média ao longo da 7ª dimensão (trials)
EVM_mean = mean(EVM, 7);
K_fixed = numel(range_snapshots);

%%  GRÁFICOS

for iRadius = 1:nRadius
    for iISR = 1:nISR
        for iSNR = 1:nSNR
            fig = figure('Name', "RMSE vs SNR | ISR=" + string(range_ISR_dB(iISR)) + " (r = " + string(range_radius(iRadius)) + ")", ...
                'NumberTitle', 'off');
            hold on; grid on;

            legend_strings = strings(0);   % inicializa como array de strings vazio

            for iCoupling = 1:nCoupling
                for m = 1:nMethods
                    plot(range_SNR_dB, squeeze(RMSE_mean(m, :, iISR, iRadius, iCoupling, K_fixed)), ...
                        markers(m), 'Color', colors(m),'LineWidth', 2, 'LineStyle', line_style(iCoupling));
                end

            end
            for m = 1:nMethods
                legend_strings(end+1) = methods(m);
            end

            xlabel("SNR (dB)",'Interpreter','latex');
            ylabel("RMSE (degrees)",'Interpreter','latex');
            yscale("log");
            ytickformat('%.2f') % 'f' forces fixed-point notation instead of scientific
            ylim([1 100]);
            % title(sprintf('RMSE vs SNR | ISR=%d dB | Snapshots=%d', ...
            %     range_ISR_dB(iISR), range_snapshots(end)));

            legend(legend_strings, 'Location', 'northwest','Interpreter','latex');
            % legend(legend_strings, 'Location', 'best','Interpreter','latex');

            fileName = sprintf(['RMSE_vs_SNR_ISR_%ddB_Radius_%' ...
                'f.png'], ...
                range_ISR_dB(iISR), range_radius(iRadius));

            exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);

            fileName = sprintf('RMSE_vs_SNR_ISR_%ddB_Radius_%f.tex', ...
                range_ISR_dB(iISR), range_radius(iRadius));

            %cleanfigure;
            % Define LaTeX macros for width and height (you will define these in your .tex file)
            matlab2tikz(fullfile(outDir, fileName), 'width', '\figurewidth', 'height', '\figureheight');
        end
    end
end


return;


%% ======= Beampatterns =======
ang_array = -180:0.01:180;    % graus

[phi_beampattern, B_dB_CAPON] = utils.beampattern_db_uca(w_capon, M, r, lambda, ang_array);

figure('Name','Beampattern','Position',[100 100 900 500]);
hold on; grid on; box on;

plot(phi_beampattern, B_dB_CAPON,'LineWidth',1.4);

xline(phi_hat_deg,'g--','phi\_hat\_deg');
% xline(phi_hat_deg_noisy,'g--','phi\_hat\_deg\_noisy');
xline(phi_hat_int_deg,'g--','phi\_hat\_int\_deg');
% xline(phi_hat_int_deg_noisy,'g--','phi\_hat\_int\_deg\_noisy');

xlabel('\phi (graus)');
ylabel('Resposta (dB)');
title('Padrão de radiação');
legend('Capon', ...
    'Location','best');

% add_test_parameters(params);

ylim([-40 5]);                 % piso visual
xlim([-180 180]);                % faixa angular total











% =======================================================================
%% Funções auxiliares
% =======================================================================


function P = uca_positions(uca)
% Retorna P = Mx3 (x,y,z) em metros

M = uca.NumElements;
R = uca.Radius;

phi0 = uca.AngleOffset;

% Tenta ser robusto: se parecer graus, converte pra rad
if abs(phi0) > 2*pi
    phi0 = deg2rad(phi0);
end

phi = phi0 + (0:M-1).' * (2*pi/M);

P = [R*cos(phi), R*sin(phi), zeros(M,1)];
end



function Ctx = compute_Ctx_for_R(fc, M, R, Z0)
% Monta UCA e calcula Z via S-parameters, depois Ctx (pág. 29)

c = 3e8;
lambda = c/fc;

% --- Elemento ---
mp = dipole;
mp.Length  = 0.5*lambda;
mp.Width  = 0.01*lambda;

% --- UCA ---
uca = circularArray;
uca.Element = mp;
uca.NumElements = M;
uca.Radius = R;

% --- S-parameters = matriz de acoplamento ---
Sobj = sparameters(uca, fc);
S_matrix = Sobj.Parameters(:,:,1);

Z_matrix = s2z(S_matrix, Z0); % conversão Z corrijida

Zg = Z0;
Zself = diag(Z_matrix);          % Z11, Z22, ..., ZNN
denom = Zself + Zg;              % (Zjj + Zg,j)

Ctx = eye(M);

% fora-diagonais: C(i,j) = Z(i,j)/(Z(j,j)+Zg(j))
for j = 1:M
    for i = 1:M
        if i ~= j
            Ctx(i,j) = Z_matrix(i,j) / denom(j);
        end
    end
end
end


   function C = build_C_circulant_uca4(c1, c2, c3, c4)
    % UCA de 8 elementos
    first_row = [1, c1, c2, c3, c4, c3, c2, c1];
    M = 8;
    C = zeros(M,M);
    for i = 1:M
        C(i,:) = circshift(first_row, [0 i-1]);
    end
end

function J = cost_mcm_circulant4(p, a, b_hat)
    c1 = p(1) + 1j*p(2);
    c2 = p(3) + 1j*p(4);
    c3 = p(5) + 1j*p(6);
    c4 = p(7) + 1j*p(8);

    C = build_C_circulant_uca4(c1, c2, c3, c4);
    b_model = C * a;

    % remove ambiguidade de ganho escalar complexo
    alpha = (b_model' * b_hat) / (b_model' * b_model);

    err = b_hat - alpha * b_model;
    J = norm(err)^2;
end