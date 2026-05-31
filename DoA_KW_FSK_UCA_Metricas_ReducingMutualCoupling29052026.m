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

% =========================================================================
% Modo de execucao: quick_run = true para testes rapidos (poucos angulos),
% false para varredura completa de producao (todos angulos x SNRs x ISRs).
% =========================================================================
quick_run = false;   % <-- alterar conforme necessidade

% =========================================================================
% Modo "limpo" para validacao do pipeline:
%   clean_simulation = true   -> zera todas as imperfeicoes (nao-circulancia,
%                                erro de sync, perturbacao do modelo).
%                                Recupera resultados teoricos: KW LS 1-dir
%                                deve dar erro Frobenius ~1e-4.
%   clean_simulation = false  -> usa os valores definidos abaixo
%                                (cenario realista com imperfeicoes).
% =========================================================================
clean_simulation = true;   % <-- alterar conforme necessidade

% =========================================================================
% Numero de ensembles (rodadas Monte Carlo por ponto de operacao).
% Para cada combinacao de parametros, geram-se n_ensemble realizacoes
% independentes (novos sinais/ruido) e as metricas sao a MEDIA das rodadas.
% n_ensemble = 1 recupera o comportamento antigo.
% =========================================================================
n_ensemble = 1;   % <-- numero de rodadas Monte Carlo por ponto

% =========================================================================
% Beamforming: qual angulo usar para apontar o beamformer (DAS/Capon).
%   false -> usa o angulo VERDADEIRO do sinal (phi_sig). Mede o limite
%            teorico da compensacao, isolando o efeito da matriz C.
%   true  -> usa o angulo ESTIMADO pelo metodo de DoA do cenario. Mede
%            o desempenho fim-a-fim realista (erro de DoA propaga p/ o BF).
% =========================================================================
use_doa_estimate_in_bf = true;   % <-- alterar conforme necessidade

% =========================================================================
% Metodo de DoA usado como REFERENCIA nos cenarios sem self-cal
% (No-MC, No-Comp, Ideal, Pert, KW-1D, KW-MD). Os cenarios de self-cal
% (SC-*) sempre usam o proprio metodo. Alterna entre 'MUSIC' e 'CAPON'.
% =========================================================================
ref_doa_method = 'CAPON';   % 'MUSIC' | 'CAPON'

% =========================================================================
% Modulacao do interferente:
%   'fm'   -> interferente FM/narrowband CE (comportamento original)
%   'fsk2' -> interferente FSK-2 co-canal (mesma modulacao do SOI; pior
%             caso, separavel apenas espacialmente)
% =========================================================================
interf_type = 'fsk2';   % 'fm' | 'fsk2'

% =========================================================================
% Controle de quais conjuntos de plots agregados gerar:
%   plot_vs_snr -> figuras com eixo X = SNR (uma por valor de ISR)
%   plot_vs_isr -> figuras com eixo X = ISR (uma por valor de SNR)
% Util para suprimir um eixo quando ha so 1 valor varrido nele.
% =========================================================================
plot_vs_snr = true;
plot_vs_isr = false;

% =========================================================================
% Passo da varredura angular do scan de DoA (graus).
% 0.1 (default antigo) -> 3601 pontos, alta precisao mas lento.
% 0.5                  -> 721 pontos, ~5x mais rapido. Erro de
%                         quantizacao do pico fica em +/-0.25 deg, muito
%                         abaixo do RMSE tipico.
% =========================================================================
doa_scan_step = 0.5;

% =========================================================================
% Modo de varredura angular:
%   'cross'      -> phi_sig e phi_int varrem o mesmo conjunto via produto
%                   cartesiano (modo classico). nPhi x nPhi simulacoes.
%   'accuracy'   -> phi_sig varre o cone fundamental [0,45) do UCA-8 e
%                   phi_int = phi_sig + 90 (separacao fixa e larga). PAREADO.
%   'resolution' -> phi_sig fixo, phi_int varre separacoes crescentes.
%                   PAREADO.
%   'aleatory'   -> phi_sig e phi_int sorteados aleatoriamente em [-180,180]
%                   a cada ponto, respeitando separacao minima. PAREADO,
%                   n_rand_angles pares por ensemble.
% =========================================================================
experiment_mode = 'aleatory';   % 'cross' | 'accuracy' | 'resolution' | 'aleatory'

% Parametros do modo 'aleatory'
n_rand_angles  = 50;    % quantos pares (phi_sig, phi_int) sortear
min_sep_deg    = 5;    % separacao angular minima entre sinal e interferente

if quick_run
    range_SNR_dB = 6;
    range_ISR_dB = -1;
    range_phi_sig = [75 100];
    range_phi_int = [];               % nao usado no modo 'cross'
    pair_mode = 'cross';
else
    range_SNR_dB = -9:3:12;
    range_ISR_dB = [-30 -6];%-6:1:-3;

    switch experiment_mode
        case 'cross'
            range_phi_sig = -180:36:180;
            range_phi_int = [];                            % nao usado
            pair_mode = 'cross';
        case 'accuracy'
            range_phi_sig = 0:5:45;                        % cone fundamental
            range_phi_int = range_phi_sig + 15;            % separacao fixa
            pair_mode = 'paired';
        case 'resolution'
            phi_sig_fixed = 30;                            % direcao "boa" fixa
            sep_grid = [5 10 15 20 30 45 60 90 120];
            range_phi_sig = repmat(phi_sig_fixed, size(sep_grid));
            range_phi_int = phi_sig_fixed + sep_grid;
            pair_mode = 'paired';
        case 'aleatory'
            % Sorteia n_rand_angles pares respeitando separacao minima.
            % Os angulos sao FIXADOS aqui (mesmos pares para todos os
            % ensembles); o que varia entre ensembles e' o sinal/ruido.
            rng(12345, 'twister');   % reprodutibilidade dos angulos
            range_phi_sig = zeros(1, n_rand_angles);
            range_phi_int = zeros(1, n_rand_angles);
            for kk = 1:n_rand_angles
                ps = -180 + 360*rand;
                pi_ = -180 + 360*rand;
                % Garante separacao minima (distancia circular)
                d = abs(ps - pi_); d = min(d, 360 - d);
                while d < min_sep_deg
                    pi_ = -180 + 360*rand;
                    d = abs(ps - pi_); d = min(d, 360 - d);
                end
                range_phi_sig(kk) = ps;
                range_phi_int(kk) = pi_;
            end
            pair_mode = 'paired';
        otherwise
            error('experiment_mode invalido: %s', experiment_mode);
    end
end

% Para compatibilidade com o resto do script (que ainda usa range_phi):
% no modo 'cross', range_phi e' o conjunto unico; nos demais, e' range_phi_sig
% (apenas para nPhi).
range_phi = range_phi_sig;
range_snapshots = 2000:3000:N_DOA;
range_radius = [0.20]; %[0.25 0.2 0.15 0.1];

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
range_coupling = [0, 1, 2, 4, 5, 6, 7, 8, 9];

% --- parametros da calibracao (Khan 2020): UMA medicao de uma direcao conhecida ---
phi_cal_deg   = 0;     % azimute da fonte de calibracao (Coupling=4)
theta_cal_deg = 90;    % elevacao (plano XY)
SNR_cal_dB    = 12;    % SNR alto (camara anecoica)
ISR_cal_dB    = -40;   % praticamente sem interferente

% --- Modelo perturbado (Coupling=3): incerteza relativa em Z_t ---
%       Simula erro de simulacao EM / variabilidade de fabricacao do PCB.
pert_level_rel = 0.05;     % 5% de incerteza relativa nos coeficientes
rng(2025, 'twister');      % reprodutibilidade

% --- KW varias direcoes (Coupling=5): P direcoes de calibracao ---
phi_cal_multi_deg = [0, 30, 60, 90];   % evita 22.5 (singular UCA-8)
P_multi = numel(phi_cal_multi_deg);

% --- IMPERFEICOES DO MUNDO REAL ---
% (a) Quebra de circulancia: matriz de acoplamento real tem desvios em
%     relacao ao modelo circulante perfeito (variabilidade de fabricacao,
%     comprimentos de cabos diferentes, etc.). Modelado como perturbacao
%     aditiva complexa.
C_noncirc_level = 0.02;        % 2% de quebra de circulancia (fixo)

% (b) Erro de sincronismo temporal: receptor correlaciona X com q assumindo
%     alinhamento perfeito; na pratica ha offset fracionario de amostras.
%     Modelado como deslocamento (shift) aplicado ao q usado pelo receptor.
sync_error_samples = 0.5;      % offset fracionario em amostras (fixo)

% --- Override: se clean_simulation=true, zera todas as imperfeicoes ---
if clean_simulation
    pert_level_rel     = 0;
    C_noncirc_level    = 0;
    sync_error_samples = 0;
    fprintf('=== MODO LIMPO ATIVADO: imperfeicoes zeradas ===\n');
    fprintf('    (pert_level_rel=0, C_noncirc_level=0, sync_error_samples=0)\n');
end

% --- Self-cal (Coupling=6..9): parametros do alternante ---
selfcal_max_iter = 12;
selfcal_grid_deg = -180:0.5:180;
selfcal_methods  = {'DAS', 'CAPON', 'MUSIC', 'KW'};   % ordem para Coupling 6..9
selfcal_damping  = 0.5;        % subrelaxacao: suaviza oscilacoes (1.0 = sem damping)

% --- Inits da self-cal: roda DUAS VEZES, uma com cada init, e duplica plots ---
range_selfcal_init = {'identity'}; %, 'kw_offline'};
nInits             = numel(range_selfcal_init);

teste_phi = 75;

methodsDoa = ["KW","DAS","MPDR","MUSIC"];
methodsBeamforming = ["SEM BF","MPDR","DAS"];
markers = ["o-","x-","s-","d-","^-","v-","*-","+-"];
colors =  ["red", "green", "blue", "black", "magenta", "cyan", "yellow"];
line_style =   ["-", "--", ":","-."]; % linha continua sem acoplamento % linha tracejada com acoplamento

% Criar pasta 'graficos'
outDir = fullfile(pwd, 'new_graficos2');
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
% C_true sera redefinido abaixo, apos aplicar a quebra de circulancia.

% Coeficientes unicos da estrutura circulante (para diagnostico e plot)
% Para M=8: c1, c2, c3 aparecem em 2 posicoes; c4 (distancia M/2) aparece em 1.
% Sao apenas referenciais; a matriz REAL do canal pode nao ser perfeitamente circulante.


Coupling_matrices = zeros(M, M, nRadius);  % matriz IDEAL circulante (referencia)
Coupling_matrices_real = zeros(M, M, nRadius); % matriz REAL aplicada no canal
% Matrizes a ser usada nas simulações
for iRadius = 1:nRadius
    radius = range_radius(iRadius)*lambda;
    % Ctx = Z;
    Ctx = compute_Ctx_for_R(fc, M, radius, Z0);
    Coupling_matrices(:,:,iRadius) = Ctx; % matriz ideal circulante

    % --- IMPERFEICAO 1: Quebra de circulancia ---
    % Adiciona perturbacao aleatoria (NAO circulante) com intensidade
    % proporcional a magnitude media dos elementos de Ctx. A matriz "real"
    % do canal e' ligeiramente nao-circulante; a estimacao continua
    % assumindo circulancia, e o desvio aparece como erro residual.
    pert_aditiva = (randn(M) + 1j*randn(M)) / sqrt(2) * ...
                   C_noncirc_level * mean(abs(Ctx(:)));
    Coupling_matrices_real(:,:,iRadius) = Ctx + pert_aditiva;
end

% Para fins de comparacao com C_hat, usamos a matriz "real" (nao-circulante)
% como ground-truth: e' essa que esta efetivamente no canal e que queremos
% inverter. A estimacao por LS circulante NUNCA conseguira recupera-la
% perfeitamente, porque assume circulancia; o erro Frobenius residual
% reflete justamente esse limite imposto pela imperfeicao do modelo.
C_true = Coupling_matrices_real(:,:,1);  % redefine para refletir realidade
c_true_vec = zeros(M/2, 1);
for k = 1:M/2
    c_true_vec(k) = Coupling_matrices(1, 1+k, 1);   % coef. ideais (p/ diagnostico)
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
    % - usa a matriz REAL (nao-circulante), pois o sinal de calibracao
    %   tambem passa pelo mesmo hardware imperfeito.
    X_cal = Coupling_matrices_real(:,:,iRadius) * Xsig_cal + ...
            Coupling_matrices_real(:,:,iRadius) * Xint_cal + Xn_cal;

    % --- 2. Estimativa da assinatura espacial via formas de onda conhecidas ---
    %  b_hat = X * q* / (q'*q)   (mesma logica do KW-DoA)
    %
    % IMPERFEICAO 2: Erro de sincronismo. O receptor "acha" que o q
    % usado na correlacao esta alinhado com o sinal recebido, mas na
    % pratica ha um offset fracionario. Modelado aplicando shift no q
    % usado pelo receptor.
    q_cal_misaligned = apply_fractional_shift(q_cal(:), sync_error_samples);
    b_hat_cal = X_cal * conj(q_cal_misaligned) / (q_cal_misaligned' * q_cal_misaligned);

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
        X_p = Coupling_matrices_real(:,:,iRadius) * Xs_p + ...
              Coupling_matrices_real(:,:,iRadius) * Xi_p + Xn_p;
        % Imperfeicao 2: q desalinhado tambem aqui
        q_p_misaligned = apply_fractional_shift(q_p(:), sync_error_samples);
        B_multi(:, pp) = X_p * conj(q_p_misaligned) / (q_p_misaligned' * q_p_misaligned);
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
% Salva resultados no caso de referencia (phi_sig=75, phi_int=100).
% MATRIZ 2D: results_cmp(iCoupling, iInit), onde iInit varia em 1..nInits
% para Coupling 6..9 (self-cal) e ignora-se iInit > 1 para Coupling 0..5.

% =========================================================================
% Arrays de agregacao estatistica (varredura completa SNR x ISR x phi_sig x phi_int).
% Indexacao: (iCoupling, iInit, iSNR, iISR, iPhi, iiPhi)
%   - iPhi  -> indice em range_phi (phi_sig)
%   - iiPhi -> indice em range_phi (phi_int)
% NOTA: com ensembles, os arrays guardam a SOMA das rodadas; a divisao
% pela contagem (media) e' feita ao final do loop principal. agg_count
% conta quantas rodadas validas escreveram em cada celula.
% =========================================================================
agg_dims = [nCoupling, nInits, nSNR, nISR, nPhi, nPhi];
agg_phi_err_KW    = zeros(agg_dims);   % SOMA de |phi_KW - phi_sig|
agg_phi_err_DAS   = zeros(agg_dims);
agg_phi_err_MVDR  = zeros(agg_dims);
agg_phi_err_MUSIC = zeros(agg_dims);
agg_eps_F         = zeros(agg_dims);   % SOMA do erro Frobenius
agg_phi_sig       = nan(agg_dims);     % phi_sig real (nao acumula; debug)
agg_phi_KW        = zeros(agg_dims);   % SOMA do phi estimado
agg_BER_DAS       = zeros(agg_dims);   % SOMA BER apos beamforming DAS
agg_BER_MVDR      = zeros(agg_dims);   % SOMA BER apos beamforming MVDR/Capon
agg_EVM_DAS       = zeros(agg_dims);   % SOMA EVM apos beamforming DAS (dB)
agg_EVM_MVDR      = zeros(agg_dims);   % SOMA EVM apos beamforming MVDR/Capon (dB)
agg_count         = zeros(agg_dims);   % nº de rodadas validas por celula

results_cmp(nCoupling, nInits) = struct( ...
    'label',[], 'init_label',[], ...
    'phi_KW',NaN, 'phi_DAS',NaN, 'phi_MVDR',NaN, 'phi_MUSIC',NaN, ...
    'P_DAS_dB',[], 'P_MVDR_dB',[], 'P_MUSIC_dB',[], 'phi_scan',[], ...
    'BP_DAS_dB',[], 'BP_Capon_dB',[], 'phi_bp',[], ...
    'BER_in',NaN, 'BER_DAS',NaN, 'BER_MVDR',NaN, ...
    'EVM_in',NaN, 'EVM_DAS',NaN, 'EVM_MVDR',NaN, ...
    'sc_history',[]);
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
    for iI = 1:nInits
        results_cmp(ic, iI).label      = labels_cmp(range_coupling(ic)+1);
        results_cmp(ic, iI).init_label = range_selfcal_init{iI};
    end
end

% =========================================================================
% Pre-computacao: matriz de steering vectors do scan de DoA, por raio.
% Como theta_sig_deg, M, lambda e o vetor phi_scan_doa sao fixos (radius
% varia por iRadius), pre-computamos A_scan_doa_cell{iRadius} fora dos
% loops, evitando milhoes de chamadas a steering_vec_uca dentro do sweep.
% =========================================================================
phi_scan_doa_global = -180:doa_scan_step:180;     % grade do DoA
n_scan_doa = numel(phi_scan_doa_global);
A_scan_doa_cell = cell(nRadius, 1);
for iRadius_pre = 1:nRadius
    radius_pre = range_radius(iRadius_pre) * lambda;
    A = zeros(M, n_scan_doa);
    for kk = 1:n_scan_doa
        A(:, kk) = utils.steering_vec_uca(M, radius_pre, lambda, ...
                                          theta_sig_deg, phi_scan_doa_global(kk));
    end
    A_scan_doa_cell{iRadius_pre} = A;
end

% =========================================================================
% Inicializa cronometro e contador de progresso
% =========================================================================
sweep_t0 = tic;
if strcmp(pair_mode, 'paired')
    TotalSim = nSNR * nISR * nRadius * nCoupling * nPhi * n_ensemble;
else
    % modo cross: descarta a diagonal phi_sig==phi_int (nPhi casos pulados)
    TotalSim = nSNR * nISR * nRadius * nCoupling * (nPhi*nPhi - nPhi) * n_ensemble;
end
iTotal   = 0;
fprintf('\n=== INICIO DO LOOP PRINCIPAL ===\n');
fprintf('Total de simulacoes: %d  (SNR x ISR x Radius x Coupling x phi_sig x phi_int)\n', TotalSim);
fprintf('================================\n\n');

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
                    phi_sig_deg = range_phi_sig(iPhi);

                    % Define range para o segundo loop conforme pair_mode:
                    %   - 'cross': phi_int varre range_phi_sig (todos)
                    %   - 'paired': phi_int e' apenas range_phi_int(iPhi)
                    if strcmp(pair_mode, 'paired')
                        iiPhi_range  = iPhi;
                        phi_int_list = range_phi_int(iPhi);
                    else
                        iiPhi_range  = 1:nPhi;
                        phi_int_list = range_phi_sig;   % varre o mesmo conjunto
                    end

                    for iiPhi_local = 1:numel(iiPhi_range)
                        iiPhi = iiPhi_range(iiPhi_local);
                        phi_int_deg = phi_int_list(iiPhi_local);

                        if(phi_int_deg ~= phi_sig_deg)
                            j=1;

                            % ===== LOOP DE ENSEMBLE (Monte Carlo) =====
                            % Para cada ponto de operacao, geram-se
                            % n_ensemble realizacoes independentes. As
                            % metricas sao acumuladas e a media e' feita
                            % ao final do loop principal.
                            for iEns = 1:n_ensemble

                            % --- Progresso da simulacao ---
                            iTotal = iTotal + 1;
                            elapsed = toc(sweep_t0);
                            if iTotal > 1
                                eta = elapsed * (TotalSim - iTotal) / (iTotal - 1);
                                eta_str = datestr(seconds(eta), 'HH:MM:SS');
                            else
                                eta_str = '--:--:--';
                            end
                            fprintf(['--> Sim %d / %d (%.1f%%)  |  ' ...
                                     'SNR=%+d  ISR=%+d  Coup=%d  r=%d  ' ...
                                     'phi_s=%+.0f  phi_i=%+.0f  ens=%d/%d  |  ' ...
                                     'decorrido %s  ETA %s\n'], ...
                                    iTotal, TotalSim, 100*iTotal/TotalSim, ...
                                    range_SNR_dB(iSNR), range_ISR_dB(iISR), ...
                                    range_coupling(iCoupling), iRadius, ...
                                    phi_sig_deg, phi_int_deg, iEns, n_ensemble, ...
                                    datestr(seconds(elapsed), 'HH:MM:SS'), eta_str);

                            %phi_sig_deg = -180 + 360*rand;
                            %phi_int_deg = -180 + 360*rand;

                            % ----- Forma de onda conhecida (q) e interferidor "ruído" (r) -----

                            [X, q, r_int, Xsig, Xint, Xn, bits, pam_rrc_tx, pam_rect, sym_tx, qn, taus_sig, taus_int] = ...
                                utils.simulate_fsk_data_uca(M, radius, lambda, phi_sig_deg, phi_int_deg, ...
                                theta_sig_deg, theta_int_deg, SNR_dB, ISR_dB, N, fs, Rs, sps, alpha, span, fd, interf_type);

                            % Aplicando o Mutual Coupling no canal
                            if Coupling == 0
                                % Caso 0: ideal, sem acoplamento
                                Coupling_matrix = eye(M);
                            else
                                % Casos 1..4: canal SEMPRE com acoplamento real (C_true).
                                % O que muda entre eles e' a estrategia de
                                % compensacao aplicada APOS a recepcao.
                                Coupling_matrix = Coupling_matrices_real(:,:,iRadius);
                            end

                            % X = Coupling_matrix * X;
                            % X = Coupling_matrix_sig * Xsig + Coupling_matrix_int * Xint + Xn;
                            X_1antenna = Xsig + Xint + Xn;
                            X_1antenna = X_1antenna(1,:);
                            Xq = Xsig + Xint + Xn;
                            % X = Coupling_matrix * Xsig + Coupling_matrix * Xint + Xn;
                            X = Coupling_matrix * Xsig + Coupling_matrix * Xint + Xn;

                            % --- Casos com compensacao (Coupling 2..9) ---
                            % O receptor compensa o acoplamento (sem saber a DoA do sinal).
                            %
                            % Loop sobre iInit: para Coupling 6..9 (self-cal),
                            % roda 2 vezes, uma para cada init (identity / kw_offline).
                            % Para Coupling 0..5, roda apenas iInit=1 (init e' irrelevante).
                            if Coupling >= 6 && Coupling <= 9
                                nInitsThis = nInits;
                            else
                                nInitsThis = 1;
                            end

                            X_baseline = X;   % salva X antes da compensacao para reusar
                            % q desalinhado: usado em todas as correlacoes do receptor
                            % (compensacao + DoA-KW operacional)
                            q_misaligned = apply_fractional_shift(q(:), sync_error_samples);
                            for iInit = 1:nInitsThis
                                X = X_baseline;   % restaura para cada init
                                sc_hist_iter = [];

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

                                    % Determina o init para esta passada
                                    init_str = range_selfcal_init{iInit};
                                    switch init_str
                                        case 'identity'
                                            init_arg = 'identity';
                                        case 'kw_offline'
                                            % Usa diretamente a matriz da
                                            % calibracao offline (Coupling 4)
                                            init_arg = C_kw_1dir{iRadius};
                                        otherwise
                                            init_arg = init_str;
                                    end

                                    % q desalinhado ja' foi calculado acima
                                    [C_sc, ~, ~, ~, sc_hist_iter] = estimate_C_selfcal(...
                                        X, q_misaligned, M, radius_m_sc, lambda, ...
                                        sc_method, selfcal_grid_deg, ...
                                        selfcal_max_iter, [], [], C_true, ...
                                        selfcal_damping, init_arg);
                                    X = (C_sc \ X);
                                end

                            K = range_snapshots(end);
                            % ----- DoA KW
                            % theta_hat_deg   = doa_kw_2007(X(:,1:K), q(1:K), M, d, lambda, 1);        % ell = 1 com M=4

                            beta = 2*pi*(0:M-1)'/M;
                            [theta_hat_deg, phi_hat_deg] = doa_kw_uca(X(:,1:K), q_misaligned(1:K).', radius, lambda, beta);
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
                            % Usa a matriz pre-computada A_scan_doa (M x n_scan_doa)
                            % e a grade global phi_scan_doa_global, evitando
                            % recalcular steering vectors a cada iteracao.
                            phi_scan = phi_scan_doa_global;
                            A_scan  = A_scan_doa_cell{iRadius};
                            theta_scan = 90;            % fixa em 90° (plano XY)

                            Rinv  = inv(Rxx_dl);
                            EnEnH = En * En';   % cache: usado em todos os angulos

                            % Vetorizacao: para A_scan (MxL) e B (MxM) hermitiana,
                            % o vetor [a_k' * B * a_k] (k=1..L) e' sum(conj(A) .* (B*A), 1).
                            % Isso elimina o loop angular e usa BLAS.
                            BA_DAS   = Rxx   * A_scan;
                            BA_MVDR  = Rinv  * A_scan;
                            BA_MUSIC = EnEnH * A_scan;
                            P_DAS   = abs( sum(conj(A_scan) .* BA_DAS,   1) );
                            denom_mvdr   = real( sum(conj(A_scan) .* BA_MVDR,  1) );
                            P_MVDR  = 1 ./ max(denom_mvdr, eps);
                            denom_music  = real( sum(conj(A_scan) .* BA_MUSIC, 1) );
                            P_MUSIC = 1 ./ max(denom_music, eps);

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

                            % --- Angulo para apontar o beamformer ---
                            % Se use_doa_estimate_in_bf=false: usa phi_sig (verdadeiro).
                            % Se true: usa o phi estimado pelo metodo do cenario.
                            if ~use_doa_estimate_in_bf
                                phi_bf_target = phi_sig_deg;
                            else
                                switch Coupling
                                    case 6, phi_bf_target = phi_DAS;    % SC-DAS
                                    case 7, phi_bf_target = phi_MVDR;   % SC-CAPON
                                    case 8, phi_bf_target = phi_MUSIC;  % SC-MUSIC
                                    case 9, phi_bf_target = phi_hat_deg;% SC-KW
                                    otherwise
                                        % Referencias: usa o metodo configurado
                                        if strcmpi(ref_doa_method, 'CAPON')
                                            phi_bf_target = phi_MVDR;
                                        else
                                            phi_bf_target = phi_MUSIC;
                                        end
                                end
                            end

                            % --- Agregacao estatistica (varredura completa) ---
                            % Calcula erro angular de cada metodo vs phi_sig real.
                            % Usa wrap [-180, 180) para evitar pulos.
                            wrap_err = @(e) mod(e + 180, 360) - 180;
                            idx6 = {iCoupling, iInit, iSNR, iISR, iPhi, iiPhi};
                            agg_phi_err_KW(idx6{:})    = agg_phi_err_KW(idx6{:})    + abs(wrap_err(phi_hat_deg - phi_sig_deg));
                            agg_phi_err_DAS(idx6{:})   = agg_phi_err_DAS(idx6{:})   + abs(wrap_err(phi_DAS - phi_sig_deg));
                            agg_phi_err_MVDR(idx6{:})  = agg_phi_err_MVDR(idx6{:})  + abs(wrap_err(phi_MVDR - phi_sig_deg));
                            agg_phi_err_MUSIC(idx6{:}) = agg_phi_err_MUSIC(idx6{:}) + abs(wrap_err(phi_MUSIC - phi_sig_deg));
                            agg_phi_sig(idx6{:})       = phi_sig_deg;   % constante; nao acumula
                            agg_phi_KW(idx6{:})        = agg_phi_KW(idx6{:})        + phi_hat_deg;
                            % conta a rodada para esta celula (por init)
                            agg_count(idx6{:}) = agg_count(idx6{:}) + 1;

                            % Erro Frobenius da matriz de acoplamento usada nesta combinacao
                            switch Coupling
                                case 2, eps_F_here = 0;   % oracle
                                case 3, eps_F_here = err_Cb_F(iRadius);
                                case 4, eps_F_here = err_C_kw1(iRadius);
                                case 5, eps_F_here = err_C_kwmd(iRadius);
                                case {6,7,8,9}
                                    if ~isempty(sc_hist_iter)
                                        eps_F_here = sc_hist_iter.err_F_per_iter(end);
                                    else
                                        eps_F_here = NaN;
                                    end
                                otherwise
                                    eps_F_here = NaN;   % Coupling 0 e 1: sem ^C
                            end
                            if ~isnan(eps_F_here)
                                agg_eps_F(idx6{:}) = agg_eps_F(idx6{:}) + eps_F_here;
                            end

                            % --- NOVO: armazenar para comparacao final
                            % (caso de referencia: phi_sig=75) ---
                            % OPTIMIZACAO: vetores pesados (P_*_dB, phi_scan)
                            % so sao usados pelos plots quick_run. No full
                            % sweep eles seriam alocados/copiados a cada
                            % ensemble sem nunca serem lidos, gerando overhead
                            % significativo. Guardamos so escalares sempre,
                            % e vetores apenas em quick_run.
                            if abs(phi_sig_deg - 75) < 1e-9
                                results_cmp(iCoupling, iInit).phi_KW    = phi_hat_deg;
                                results_cmp(iCoupling, iInit).phi_DAS   = phi_DAS;
                                results_cmp(iCoupling, iInit).phi_MVDR  = phi_MVDR;
                                results_cmp(iCoupling, iInit).phi_MUSIC = phi_MUSIC;
                                if quick_run
                                    results_cmp(iCoupling, iInit).P_DAS_dB   = P_DAS_dB;
                                    results_cmp(iCoupling, iInit).P_MVDR_dB  = P_MVDR_dB;
                                    results_cmp(iCoupling, iInit).P_MUSIC_dB = P_MUSIC_dB;
                                    results_cmp(iCoupling, iInit).phi_scan   = phi_scan;
                                end
                                % Para self-cal (Coupling 6..9), salva history
                                if ~isempty(sc_hist_iter)
                                    results_cmp(iCoupling, iInit).sc_history = sc_hist_iter;
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
                            % Garante que o angulo-alvo do BF existe EXATAMENTE
                            % na grade (pode ser phi_sig ou um phi estimado
                            % fracionario). Sem isso, a condicao de igualdade
                            % na gravacao do BER/EVM nunca dispararia.
                            if ~ismember(phi_bf_target, phi_scan)
                                phi_scan = sort([phi_scan, phi_bf_target]);
                            end
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
                                    results_cmp(iCoupling, iInit).phi_bp      = phi_beampattern;
                                    results_cmp(iCoupling, iInit).BP_DAS_dB   = B_dB_DAS;
                                    results_cmp(iCoupling, iInit).BP_Capon_dB = B_dB_Capon;
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

                                % fprintf('BER: RX %.2f%% / Capon %.2f%% / DAS %.2f%% \n', ...
                                %     BER_in*100, BER_mvdr*100, BER_das*100);

                                % Calcula EVM
                                [EVM_in,  EVMdB_in]                 = utils.calc_evm_real(sym_rx1,  sym_tx);
                                [EVM_das, EVMdB_das]                = utils.calc_evm_real(sym_rx2, sym_tx);
                                [EVM_mvdr,EVMdB_mvdr]               = utils.calc_evm_real(sym_rx3, sym_tx);

                                EVM_scan_in(ang) = 20*log10(EVM_in);
                                EVM_scan_mvdr(ang) = 20*log10(EVM_mvdr);
                                EVM_scan_das(ang) = 20*log10(EVM_das);
                                 
                                % fprintf('EVM (dB): Rx: %.2f | Capon: %.2f | DAS: %.2f\n', ...
                                %     EVMdB_in, EVMdB_mvdr, EVMdB_das);

                                % --- NOVO: salvar BER/EVM para comparacao final entre cenarios ---
                                if abs(phi_scan(ang) - teste_phi) < 1e-9 && abs(phi_sig_deg - 75) < 1e-9
                                    results_cmp(iCoupling, iInit).BER_in   = BER_in   * 100;
                                    results_cmp(iCoupling, iInit).BER_DAS  = BER_das  * 100;
                                    results_cmp(iCoupling, iInit).BER_MVDR = BER_mvdr * 100;
                                    results_cmp(iCoupling, iInit).EVM_in   = 20*log10(EVM_in);
                                    results_cmp(iCoupling, iInit).EVM_DAS  = 20*log10(EVM_das);
                                    results_cmp(iCoupling, iInit).EVM_MVDR = 20*log10(EVM_mvdr);
                                end

                                % --- Agregacao de BER/EVM para todas as combinacoes ---
                                % Salva quando o beamformer aponta para o angulo-alvo
                                % (phi_sig verdadeiro OU phi estimado, conforme a flag
                                %  use_doa_estimate_in_bf).
                                if abs(phi_scan(ang) - phi_bf_target) < 1e-9
                                    idxB = {iCoupling, iInit, iSNR, iISR, iPhi, iiPhi};
                                    agg_BER_DAS(idxB{:})  = agg_BER_DAS(idxB{:})  + BER_das  * 100;
                                    agg_BER_MVDR(idxB{:}) = agg_BER_MVDR(idxB{:}) + BER_mvdr * 100;
                                    agg_EVM_DAS(idxB{:})  = agg_EVM_DAS(idxB{:})  + 20*log10(EVM_das);
                                    agg_EVM_MVDR(idxB{:}) = agg_EVM_MVDR(idxB{:}) + 20*log10(EVM_mvdr);
                                end

                            end

                            end   % --- fim do for iInit ---

                            end   % --- fim do for iEns (ensemble) ---


                            if abs(phi_sig_deg - 75) < 1e-9 && quick_run
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

%%
fprintf('\n=== LOOP PRINCIPAL CONCLUIDO em %s ===\n', ...
        datestr(seconds(toc(sweep_t0)), 'HH:MM:SS'));

% =========================================================================
% MEDIA DOS ENSEMBLES: divide a soma acumulada pelo numero de rodadas.
% Celulas com contagem zero (combinacoes nao executadas, ex: phi_sig==phi_int)
% viram NaN para serem ignoradas na agregacao posterior ('omitnan').
% =========================================================================
cnt = agg_count;
cnt(cnt == 0) = NaN;   % evita divisao por zero -> NaN
agg_phi_err_KW    = agg_phi_err_KW    ./ cnt;
agg_phi_err_DAS   = agg_phi_err_DAS   ./ cnt;
agg_phi_err_MVDR  = agg_phi_err_MVDR  ./ cnt;
agg_phi_err_MUSIC = agg_phi_err_MUSIC ./ cnt;
agg_phi_KW        = agg_phi_KW        ./ cnt;
agg_eps_F         = agg_eps_F         ./ cnt;
agg_BER_DAS       = agg_BER_DAS       ./ cnt;
agg_BER_MVDR      = agg_BER_MVDR      ./ cnt;
agg_EVM_DAS       = agg_EVM_DAS       ./ cnt;
agg_EVM_MVDR      = agg_EVM_MVDR      ./ cnt;
fprintf('Media de %d ensemble(s) aplicada.\n', n_ensemble);
fprintf('Gerando plots...\n\n');

% =========================================================================
% PLOT FINAL: Comparacao dos cenarios (sem acoplamento / com / compensado)
% phi_sig = 75 deg, phi_int = 100 deg, teste_phi = 75
% =========================================================================
phi_sig_ref = 75;
phi_int_ref = 100;
teste_phi_ref = teste_phi;

% --- Replica resultados dos cenarios sem dependencia de init ---
% (Coupling 0..5 sao identicos para qualquer iInit; so foram populados em iInit=1.
%  Para o plot ficar consistente, copia para os demais iInits.)
for ic = 1:nCoupling
    coup_code = range_coupling(ic);
    if coup_code < 6 || coup_code > 9
        for iI = 2:nInits
            results_cmp(ic, iI) = results_cmp(ic, 1);
            results_cmp(ic, iI).init_label = range_selfcal_init{iI};
        end
    end
end

% --- Identifica quais cenarios foram efetivamente preenchidos (por init) ---
valid = false(nCoupling, nInits);
for ic = 1:nCoupling
    for iI = 1:nInits
        valid(ic, iI) = ~isempty(results_cmp(ic, iI).phi_scan) && ...
                       ~isnan(results_cmp(ic, iI).phi_KW);
    end
end

% Loop sobre os inits: gera um conjunto completo de plots por init
% (somente em quick_run; no full sweep esses plots por cenario nao sao uteis)
if quick_run
for iInitPlot = 1:nInits
    init_tag = range_selfcal_init{iInitPlot};
    ic_valid = find(valid(:, iInitPlot));

    if ~isempty(ic_valid)
        % --- Tabela de estimativas de DoA com erro absoluto ---
        fprintf('\n=========================================================\n');
        fprintf('Comparacao final entre cenarios (init=%s, phi_sig=%g, phi_int=%g)\n', ...
                init_tag, phi_sig_ref, phi_int_ref);
        fprintf('=========================================================\n');
    fprintf('%-25s | %8s | %8s | %8s | %8s\n', ...
            'Cenario','KW','DAS','MVDR','MUSIC');
    fprintf('---------------------------------------------------------\n');
    for k = 1:numel(ic_valid)
        ic = ic_valid(k);
        fprintf('%-25s | %8.3f | %8.3f | %8.3f | %8.3f\n', ...
            char(results_cmp(ic, iInitPlot).label), ...
            results_cmp(ic, iInitPlot).phi_KW, results_cmp(ic, iInitPlot).phi_DAS, ...
            results_cmp(ic, iInitPlot).phi_MVDR, results_cmp(ic, iInitPlot).phi_MUSIC);
    end
    fprintf('---------------------------------------------------------\n');
    fprintf('Erro absoluto (graus) vs phi_sig = %g:\n', phi_sig_ref);
    fprintf('%-25s | %8s | %8s | %8s | %8s\n', ...
            'Cenario','|eKW|','|eDAS|','|eMVDR|','|eMUSIC|');
    fprintf('---------------------------------------------------------\n');
    for k = 1:numel(ic_valid)
        ic = ic_valid(k);
        fprintf('%-25s | %8.3f | %8.3f | %8.3f | %8.3f\n', ...
            char(results_cmp(ic, iInitPlot).label), ...
            abs(results_cmp(ic, iInitPlot).phi_KW    - phi_sig_ref), ...
            abs(results_cmp(ic, iInitPlot).phi_DAS   - phi_sig_ref), ...
            abs(results_cmp(ic, iInitPlot).phi_MVDR  - phi_sig_ref), ...
            abs(results_cmp(ic, iInitPlot).phi_MUSIC - phi_sig_ref));
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
        key = char(results_cmp(ic, iInitPlot).label);
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
        plot(results_cmp(ic, iInitPlot).phi_scan, results_cmp(ic, iInitPlot).P_MUSIC_dB, ...
             'LineWidth', 1.6, 'Color', cmap(k,:), ...
             'DisplayName', char(results_cmp(ic, iInitPlot).label));
    end
    xline(phi_sig_ref, 'k--', 'phi_{sig}', 'LabelVerticalAlignment','bottom');
    xline(phi_int_ref, 'r--', 'phi_{int}', 'LabelVerticalAlignment','bottom');
    xlabel('\phi (graus)'); ylabel('Pseudoespectro MUSIC (dB)');
    title('MUSIC: comparacao entre cenarios');
    legend('Location','best'); xlim([min(results_cmp(ic_valid(1), iInitPlot).phi_scan) ...
                                     max(results_cmp(ic_valid(1), iInitPlot).phi_scan)]);

    % (2) Beampatterns Capon sobrepostos
    subplot(2,2,2); hold on; grid on;
    for k = 1:numel(ic_valid)
        ic = ic_valid(k);
        if ~isempty(results_cmp(ic, iInitPlot).BP_Capon_dB)
            plot(results_cmp(ic, iInitPlot).phi_bp, results_cmp(ic, iInitPlot).BP_Capon_dB, ...
                 'LineWidth', 1.6, 'Color', cmap(k,:), ...
                 'DisplayName', char(results_cmp(ic, iInitPlot).label));
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
        BER_mat(k,:) = [results_cmp(ic, iInitPlot).BER_in, ...
                        results_cmp(ic, iInitPlot).BER_DAS, ...
                        results_cmp(ic, iInitPlot).BER_MVDR];
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
        EVM_mat(k,:) = [results_cmp(ic, iInitPlot).EVM_in, ...
                        results_cmp(ic, iInitPlot).EVM_DAS, ...
                        results_cmp(ic, iInitPlot).EVM_MVDR];
    end
    bh2 = bar(EVM_mat, 'grouped');
    set(gca, 'XTickLabel', short_labels, ...
        'XTickLabelRotation', 25);
    grid on;
    ylabel('EVM (dB)');
    title('EVM por cenario × beamformer');
    legend({'Rx (1 antena)','DAS','MVDR/Capon'}, 'Location','best');

    sgtitle(sprintf(['Comparacao: 10 cenarios de acoplamento (init self-cal: %s)\n' ...
                     '\\phi_{sig}=%g°, \\phi_{int}=%g°, SNR=%g dB, ISR=%g dB, ' ...
                     'raio=%.2f \\lambda, pert. modelo=%.0f%%, KW Multi P=%d'], ...
            init_tag, phi_sig_ref, phi_int_ref, range_SNR_dB(end), range_ISR_dB(end), ...
            range_radius(end), 100*pert_level_rel, P_multi));

    exportgraphics(fig_cmp, fullfile(outDir, ...
        sprintf('comparacao_cenarios_acoplamento_init_%s.png', init_tag)), ...
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
        if coup_code >= 6 && coup_code <= 9 && ~isempty(results_cmp(ic, iInitPlot).sc_history)
            sc_indices_in_results(end+1) = ic; %#ok<AGROW>
            sc_method_names{end+1}       = selfcal_methods{coup_code - 5}; %#ok<AGROW>
        end
    end

    if ~isempty(sc_indices_in_results)
        cmap_sc = lines(numel(sc_indices_in_results));

        % =====================================================================
        % FIGURA 1: Curva de erro Frobenius da self-cal
        % =====================================================================
        fig_frob = figure('Name','Self-cal: erro Frobenius vs iteracao', ...
                          'NumberTitle','off', 'Position',[100 100 1200 700]);

        % --- (1) Erro Frobenius vs iteracao -------------------------------
        subplot(1,2,1); hold on; grid on;
        for kk = 1:numel(sc_indices_in_results)
            ic = sc_indices_in_results(kk);
            h = results_cmp(ic, iInitPlot).sc_history;
            it_axis = 1:numel(h.err_F_per_iter);
            plot(it_axis, h.err_F_per_iter, 'o-', ...
                 'LineWidth', 1.8, 'Color', cmap_sc(kk,:), ...
                 'MarkerFaceColor', cmap_sc(kk,:), ...
                 'DisplayName', sprintf('SC-%s', sc_method_names{kk}));
        end
        % Linhas de referencia: erros offline
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
        title('Erro Frobenius vs iteracao');
        legend('Location','northeast');

        % --- (2) |C_t - C_{t-1}| (criterio interno de convergencia) -------
        subplot(1,2,2); hold on; grid on;
        for kk = 1:numel(sc_indices_in_results)
            ic = sc_indices_in_results(kk);
            h = results_cmp(ic, iInitPlot).sc_history;
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

        sgtitle(sprintf("Self-cal (init=%s): convergencia da matriz de acoplamento  |  \\phi_{sig}=%g°, \\phi_{int}=%g°", init_tag, phi_sig_ref, phi_int_ref), "FontWeight","bold");

        exportgraphics(fig_frob, fullfile(outDir, ...
            sprintf('selfcal_frobenius_init_%s.png', init_tag)), ...
                       'Resolution', 200);

        % =====================================================================
        % FIGURA 2: Curva de DoA estimada da self-cal
        % =====================================================================
        fig_doa = figure('Name','Self-cal: DoA estimada vs iteracao', ...
                         'NumberTitle','off', 'Position',[100 100 900 600]);

        hold on; grid on;
        for kk = 1:numel(sc_indices_in_results)
            ic = sc_indices_in_results(kk);
            h = results_cmp(ic, iInitPlot).sc_history;
            it_axis = 1:numel(h.phi_per_iter);
            plot(it_axis, h.phi_per_iter, 'o-', ...
                 'LineWidth', 1.6, 'Color', cmap_sc(kk,:), ...
                 'MarkerFaceColor', cmap_sc(kk,:), ...
                 'DisplayName', sprintf('SC-%s', sc_method_names{kk}));
        end
        yline(phi_sig_ref, 'k--', 'LineWidth', 1.4, 'DisplayName','\phi_{sig}');
        xlabel('iteracao'); ylabel('\phi estimada (°)');
        title(sprintf('Self-cal (init=%s): DoA estimada por iteracao  |  \\phi_{sig}=%g°, \\phi_{int}=%g°', ...
              init_tag, phi_sig_ref, phi_int_ref));
        legend('Location','best');

        % --- Imprime tabela resumo no terminal (em vez de no plot) --------
        fprintf('\n--- Resumo Self-Cal (init=%s) ---\n', init_tag);
        fprintf('  %-10s | %5s | %s\n', 'metodo', 'n_it', 'err.Frob.final');
        fprintf('  %s\n', repmat('-', 1, 45));
        for kk = 1:numel(sc_indices_in_results)
            ic = sc_indices_in_results(kk);
            h  = results_cmp(ic, iInitPlot).sc_history;
            fprintf('  %-10s | %5d | %.3e\n', ...
                sprintf('SC-%s', sc_method_names{kk}), h.n_iter, ...
                h.err_F_per_iter(end));
        end
        fprintf('  Referencias offline:\n');
        fprintf('    KW LS 1-dir : %.3e\n', err_C_kw1(1));
        fprintf('    KW LS Multi : %.3e\n', err_C_kwmd(1));
        fprintf('    Modelo pert.: %.3e\n\n', err_Cb_F(1));

        exportgraphics(fig_doa, fullfile(outDir, ...
            sprintf('selfcal_doa_init_%s.png', init_tag)), ...
                       'Resolution', 200);

        % =====================================================================
        % FIGURA 3: Scatter "DoA error vs Frobenius" para todos os cenarios
        % Cruza explicitamente os dois eixos: |phi_hat - phi_sig| vs erro
        % Frobenius. Posicao no plano caracteriza qualidade do cenario:
        %   - Canto inferior esquerdo: ambos bons (ideal)
        %   - Canto inferior direito : DoA boa, C ruim (PERIGOSO p/ beamforming)
        %   - Canto superior direito : ambos ruins
        % =====================================================================
        fig_scatter = figure('Name', 'DoA error vs Frobenius (todos os cenarios)', ...
                             'NumberTitle','off', 'Position',[100 100 900 700]);
        hold on; grid on;

        % Coleta pontos: percorre todos os cenarios validos (Coupling 2..9)
        scatter_points = [];   % [doa_err, frob_err, ic, label_idx]
        for ic_s = 1:nCoupling
            coup_code = range_coupling(ic_s);
            if coup_code < 2, continue; end   % pula No-MC e No-Comp (sem ^C)
            r = results_cmp(ic_s, iInitPlot);
            if isnan(r.phi_KW), continue; end

            % Erro de DoA: usa o estimador KW como referencia (mais robusto)
            doa_err_kw = abs(r.phi_KW - phi_sig_ref);
            doa_err_das = abs(r.phi_DAS - phi_sig_ref);
            doa_err_mvdr = abs(r.phi_MVDR - phi_sig_ref);

            % Erro Frobenius: para cenarios offline (2..5), usa o err calculado
            % na fase de calibracao; para self-cal (6..9), usa o ultimo do history.
            if coup_code == 2
                frob_err = 0;   % oracle por definicao
            elseif coup_code == 3
                frob_err = err_Cb_F(1);
            elseif coup_code == 4
                frob_err = err_C_kw1(1);
            elseif coup_code == 5
                frob_err = err_C_kwmd(1);
            elseif coup_code >= 6 && coup_code <= 9
                if isempty(r.sc_history), continue; end
                frob_err = r.sc_history.err_F_per_iter(end);
            end

            % Usa o KW como DoA representativa por ser mais estavel
            scatter_points(end+1, :) = [doa_err_kw, frob_err, coup_code]; %#ok<AGROW>
        end

        % Plot por categoria com cores e marcadores distintos
        markers = {'o', 's', 'd', '^', 'v', '<', '>', 'p'};
        legend_handles = [];
        legend_labels  = {};
        for kk = 1:size(scatter_points, 1)
            cc = scatter_points(kk, 3);
            switch cc
                case 2, mk = 'p'; col = [0.0 0.6 0.0]; lbl = 'Modelo ideal';
                case 3, mk = 'h'; col = [0.5 0.5 0.0]; lbl = 'Modelo perturb.';
                case 4, mk = 'd'; col = [0.0 0.4 0.7]; lbl = 'KW LS 1-dir';
                case 5, mk = 's'; col = [0.0 0.7 0.7]; lbl = 'KW LS Multi';
                case 6, mk = 'o'; col = [0.0 0.4 0.7]; lbl = 'SC-DAS';
                case 7, mk = 'o'; col = [0.85 0.3 0.1]; lbl = 'SC-CAPON';
                case 8, mk = 'o'; col = [0.93 0.69 0.13]; lbl = 'SC-MUSIC';
                case 9, mk = 'o'; col = [0.49 0.18 0.56]; lbl = 'SC-KW';
            end
            h = scatter(scatter_points(kk,1), max(scatter_points(kk,2), 1e-6), ...
                       150, mk, 'filled', 'MarkerFaceColor', col, ...
                       'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
                       'DisplayName', lbl);
            legend_handles(end+1) = h; %#ok<AGROW>
            legend_labels{end+1}  = lbl; %#ok<AGROW>
        end

        % Linhas de referencia
        xline(1, 'k:', '1° de erro DoA');
        yline(0.1, 'k:', '\epsilon_F = 0.1');
        yline(1, 'r--', '\epsilon_F = 1 (sem comp.)');

        set(gca, 'YScale', 'log');
        xlabel('|\phi_{KW} - \phi_{sig}|  (graus)');
        ylabel('||C_{hat} - C_{true}||_F / ||C_{true}||_F');
        title(sprintf('Diagnostico cruzado: erro de DoA vs erro de C (init=%s)', init_tag));
        legend(legend_handles, legend_labels, 'Location', 'best');
        xlim([-0.5 max(20, max(scatter_points(:,1))+1)]);
        exportgraphics(fig_scatter, fullfile(outDir, ...
            sprintf('scatter_doa_vs_frobenius_init_%s.png', init_tag)), ...
                       'Resolution', 200);

        % =====================================================================
        % FIGURA 4: BER/EVM apenas dos self-cals -- mostra que C ruim
        % degrada beamforming mesmo com DoA boa.
        % =====================================================================
        % Coleta apenas Coupling 6..9 para foco
        sc_for_ber = [];   % [ic, das_ber, mvdr_ber, doa_err]
        for kk = 1:numel(sc_indices_in_results)
            ic = sc_indices_in_results(kk);
            r = results_cmp(ic, iInitPlot);
            if isnan(r.BER_DAS), continue; end
            sc_for_ber(end+1, :) = [ic, r.BER_DAS, r.BER_MVDR, ...
                                    r.EVM_DAS, r.EVM_MVDR, ...
                                    abs(r.phi_KW - phi_sig_ref), ...
                                    r.sc_history.err_F_per_iter(end)]; %#ok<AGROW>
        end

        if ~isempty(sc_for_ber)
            fig_ber_sc = figure('Name', 'Self-cal: BER/EVM vs erro de C', ...
                                'NumberTitle','off', 'Position',[100 100 1300 600]);

            % --- (1) BER vs erro Frobenius ---
            subplot(1,2,1); hold on; grid on;
            for kk = 1:size(sc_for_ber, 1)
                ic = sc_for_ber(kk, 1);
                method = sc_method_names{find(sc_indices_in_results==ic,1)};
                col_idx = find(sc_indices_in_results==ic, 1);
                col = cmap_sc(col_idx, :);
                % BER do MVDR/Capon (mais sensivel a C errado)
                semilogx(max(sc_for_ber(kk,7), 1e-6), sc_for_ber(kk, 3), ...
                    'o', 'MarkerSize', 12, 'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
                    'DisplayName', sprintf('SC-%s (MVDR)', method));
                semilogx(max(sc_for_ber(kk,7), 1e-6), sc_for_ber(kk, 2), ...
                    's', 'MarkerSize', 12, 'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
                    'DisplayName', sprintf('SC-%s (DAS)', method));
            end
            set(gca, 'XScale', 'log');
            xlabel('||C_{hat} - C_{true}||_F / ||C_{true}||_F');
            ylabel('BER (%)');
            title('BER apos beamforming vs erro de C');
            xline(0.1, 'k:', '\epsilon_F = 0.1');
            xline(1, 'r--', '\epsilon_F = 1');
            legend('Location','best');

            % --- (2) EVM vs erro Frobenius ---
            subplot(1,2,2); hold on; grid on;
            for kk = 1:size(sc_for_ber, 1)
                ic = sc_for_ber(kk, 1);
                method = sc_method_names{find(sc_indices_in_results==ic,1)};
                col_idx = find(sc_indices_in_results==ic, 1);
                col = cmap_sc(col_idx, :);
                semilogx(max(sc_for_ber(kk,7), 1e-6), sc_for_ber(kk, 5), ...
                    'o', 'MarkerSize', 12, 'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
                    'DisplayName', sprintf('SC-%s (MVDR)', method));
                semilogx(max(sc_for_ber(kk,7), 1e-6), sc_for_ber(kk, 4), ...
                    's', 'MarkerSize', 12, 'MarkerFaceColor', col, ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
                    'DisplayName', sprintf('SC-%s (DAS)', method));
            end
            set(gca, 'XScale', 'log');
            xlabel('||C_{hat} - C_{true}||_F / ||C_{true}||_F');
            ylabel('EVM (dB)');
            title('EVM apos beamforming vs erro de C');
            xline(0.1, 'k:', '\epsilon_F = 0.1');
            xline(1, 'r--', '\epsilon_F = 1');
            legend('Location','best');

            sgtitle(sprintf("Self-cal (init=%s): impacto do erro de C no beamforming (circulo=MVDR, quadrado=DAS)", init_tag), "FontWeight","bold");

            exportgraphics(fig_ber_sc, fullfile(outDir, ...
                sprintf('selfcal_beamforming_vs_C_init_%s.png', init_tag)), ...
                           'Resolution', 200);
        end
    end   % --- fim do if ~isempty(sc_indices_in_results) ---

    else
        warning('Nenhum cenario foi preenchido em results_cmp para init=%s. Verifique se phi_sig=75 e teste_phi=75 estao no range simulado.', ...
                init_tag);
    end   % --- fim do if ~isempty(ic_valid) ---
end   % --- fim do for iInitPlot ---
end   % --- fim do if quick_run ---

% =========================================================================
% PLOTS AGREGADOS (so faz sentido em modo full sweep, mas funciona em ambos)
% Para cada metrica, gera duas figuras (vs SNR e vs ISR), por iInit.
% =========================================================================

% Cores e estilos por cenario (10 cenarios)
cmap_cen = [
    0.20 0.20 0.20;   % 0 No-MC      preto
    0.85 0.20 0.20;   % 1 No-Comp    vermelho
    0.10 0.65 0.10;   % 2 Ideal      verde
    0.60 0.45 0.10;   % 3 Pert       marrom-mostarda
    0.10 0.30 0.85;   % 4 KW-1D      azul
    0.20 0.75 0.80;   % 5 KW-MD      ciano
    0.95 0.50 0.10;   % 6 SC-DAS     laranja
    0.55 0.10 0.65;   % 7 SC-CAPON   roxo
    1.00 0.85 0.15;   % 8 SC-MUSIC   amarelo claro
    0    0    0;   % 9 SC-KW      rosa //indigo escuro (bem distinto do roxo) 0.30 0.00 0.50;
];
styles_cen = {':',':','--','--','-.','-.','-','-','-','-'};
markers_cen = {'none','none','s','d','d','s','o','o','o','o'};
short_labels_cen = {'No-MC', 'No-Comp', 'Ideal', 'Pert', ...
                    'KW-1D', 'KW-MD', 'SC-DAS', 'SC-CAPON', 'SC-MUSIC', 'SC-KW'};

% Lista de metodos de DoA a serem agregados
doa_methods = {'KW', 'DAS', 'MVDR', 'MUSIC'};
agg_doa = {agg_phi_err_KW, agg_phi_err_DAS, agg_phi_err_MVDR, agg_phi_err_MUSIC};

% --- Replica arrays agregados de cenarios offline (Coupling 2..5) para
%     iInit > 1, pois eles independem do init (sao calibrados antes). ---
for ic_rep = 1:nCoupling
    coup_code_rep = range_coupling(ic_rep);
    if coup_code_rep >= 2 && coup_code_rep <= 5
        for iI = 2:nInits
            agg_eps_F(ic_rep, iI, :, :, :, :)         = agg_eps_F(ic_rep, 1, :, :, :, :);
            agg_phi_err_KW(ic_rep, iI, :, :, :, :)    = agg_phi_err_KW(ic_rep, 1, :, :, :, :);
            agg_phi_err_DAS(ic_rep, iI, :, :, :, :)   = agg_phi_err_DAS(ic_rep, 1, :, :, :, :);
            agg_phi_err_MVDR(ic_rep, iI, :, :, :, :)  = agg_phi_err_MVDR(ic_rep, 1, :, :, :, :);
            agg_phi_err_MUSIC(ic_rep, iI, :, :, :, :) = agg_phi_err_MUSIC(ic_rep, 1, :, :, :, :);
            agg_BER_DAS(ic_rep, iI, :, :, :, :)       = agg_BER_DAS(ic_rep, 1, :, :, :, :);
            agg_BER_MVDR(ic_rep, iI, :, :, :, :)      = agg_BER_MVDR(ic_rep, 1, :, :, :, :);
            agg_EVM_DAS(ic_rep, iI, :, :, :, :)       = agg_EVM_DAS(ic_rep, 1, :, :, :, :);
            agg_EVM_MVDR(ic_rep, iI, :, :, :, :)      = agg_EVM_MVDR(ic_rep, 1, :, :, :, :);
        end
    end
end

% Piso visual: substitui zero/NaN por este valor no plot Frobenius em escala log
% (cenario Oracle tem eps_F = 0 que sumiria do grafico).
EPS_FLOOR = 1e-3;

% Funcoes de agregacao: RMSE sobre dimensoes especificadas
% (descarta NaN automaticamente; ex.: combinacoes phi_sig == phi_int)
rmse_along = @(A, dims) sqrt(squeeze(mean(A.^2, dims, 'omitnan')));

% --- Loop sobre inits para gerar conjunto completo ---
for iInitPlot_agg = 1:nInits
    init_tag = range_selfcal_init{iInitPlot_agg};

    % --- Mapeamento "cenario -> metodo de DoA representativo" ---
    % Self-cal (6..9): usa o proprio metodo (KW, DAS, MVDR, MUSIC)
    % Demais cenarios (0..5): nao tem metodo de DoA inerente, usam MUSIC
    %   como referencia comum.
    % Indice do metodo de referencia conforme flag ref_doa_method
    if strcmpi(ref_doa_method, 'CAPON')
        ref_midx = 3;   % MVDR/Capon
    else
        ref_midx = 4;   % MUSIC
    end
    method_idx = zeros(nCoupling, 1);   % 1=KW, 2=DAS, 3=MVDR, 4=MUSIC
    for ic_m = 1:nCoupling
        cc = range_coupling(ic_m);
        switch cc
            case 6, method_idx(ic_m) = 2;   % SC-DAS    usa DAS
            case 7, method_idx(ic_m) = 3;   % SC-CAPON  usa MVDR
            case 8, method_idx(ic_m) = 4;   % SC-MUSIC  usa MUSIC
            case 9, method_idx(ic_m) = 1;   % SC-KW     usa KW
            otherwise, method_idx(ic_m) = ref_midx;   % refs: MUSIC ou CAPON
        end
    end

    % --- Pre-calcula matrizes (SNR x ISR) por cenario, agregando sobre angulos ---
    rmse_doa_per_cen = cell(nCoupling, 1);
    rmse_eps_per_cen = cell(nCoupling, 1);
    for ic_m = 1:nCoupling
        switch method_idx(ic_m)
            case 1, A_doa = agg_phi_err_KW;
            case 2, A_doa = agg_phi_err_DAS;
            case 3, A_doa = agg_phi_err_MVDR;
            case 4, A_doa = agg_phi_err_MUSIC;
        end
        % A_doa: (Coup, Init, SNR, ISR, PhiSig, PhiInt). Agrega sobre dims 5,6.
        % Evita squeeze: extrai o slab 4D (1,1,SNR,ISR,PhiSig,PhiInt) e usa
        % mean(., [5 6]) ANTES de reshape, garantindo dimensoes consistentes
        % mesmo quando nSNR=1 ou nISR=1.
        slab_doa = A_doa(ic_m, iInitPlot_agg, :, :, :, :);
        m_doa = mean(slab_doa.^2, [5 6], 'omitnan');   % (1,1,SNR,ISR,1,1)
        rmse_doa_per_cen{ic_m} = sqrt(reshape(m_doa, [nSNR, nISR]));

        slab_eps = agg_eps_F(ic_m, iInitPlot_agg, :, :, :, :);
        m_eps = mean(slab_eps.^2, [5 6], 'omitnan');
        rmse_eps_per_cen{ic_m} = sqrt(reshape(m_eps, [nSNR, nISR]));
    end

    % =====================================================================
    % RMSE de DoA: vs SNR (uma figura por ISR), vs ISR (uma por SNR)
    % =====================================================================

    for iISR_fix = 1:nISR
        fig_h = figure('Name', sprintf('RMSE DoA vs SNR (ISR=%+d dB, init=%s)', ...
                       range_ISR_dB(iISR_fix), init_tag), ...
                       'NumberTitle','off', 'Position',[100 100 900 600]);
        hold on; grid on;
        for ic = 1:nCoupling
            coup_code = range_coupling(ic);
            curve = rmse_doa_per_cen{ic}(:, iISR_fix);
            method_label = doa_methods{method_idx(ic)};
            plot(range_SNR_dB, curve, ...
                'LineStyle', styles_cen{coup_code+1}, ...
                'Marker', markers_cen{coup_code+1}, ...
                'Color', cmap_cen(coup_code+1, :), ...
                'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'DisplayName', sprintf('%s (%s)', short_labels_cen{coup_code+1}, method_label));
        end
        set(gca, 'YScale', 'log');
        xlabel('SNR (dB)'); ylabel('RMSE de \phi (graus)');
        title(sprintf('RMSE de DoA vs SNR  |  ISR = %+d dB  |  init self-cal: %s', ...
              range_ISR_dB(iISR_fix), init_tag));
        legend('Location','best','NumColumns',2);

        exportgraphics(fig_h, fullfile(outDir, ...
            sprintf('rmse_doa_vs_SNR_ISR_%+d_init_%s.png', range_ISR_dB(iISR_fix), init_tag)), ...
            'Resolution', 200);
    end

    if plot_vs_isr
    for iSNR_fix = 1:nSNR
        fig_h = figure('Name', sprintf('RMSE DoA vs ISR (SNR=%+d dB, init=%s)', ...
                       range_SNR_dB(iSNR_fix), init_tag), ...
                       'NumberTitle','off', 'Position',[100 100 900 600]);
        hold on; grid on;
        for ic = 1:nCoupling
            coup_code = range_coupling(ic);
            curve = rmse_doa_per_cen{ic}(iSNR_fix, :);
            method_label = doa_methods{method_idx(ic)};
            plot(range_ISR_dB, curve, ...
                'LineStyle', styles_cen{coup_code+1}, ...
                'Marker', markers_cen{coup_code+1}, ...
                'Color', cmap_cen(coup_code+1, :), ...
                'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'DisplayName', sprintf('%s (%s)', short_labels_cen{coup_code+1}, method_label));
        end
        set(gca, 'YScale', 'log');
        xlabel('ISR (dB)'); ylabel('RMSE de \phi (graus)');
        title(sprintf('RMSE de DoA vs ISR  |  SNR = %+d dB  |  init self-cal: %s', ...
              range_SNR_dB(iSNR_fix), init_tag));
        legend('Location','best','NumColumns',2);

        exportgraphics(fig_h, fullfile(outDir, ...
            sprintf('rmse_doa_vs_ISR_SNR_%+d_init_%s.png', range_SNR_dB(iSNR_fix), init_tag)), ...
            'Resolution', 200);
    end
    end   % --- if plot_vs_isr ---

    % =====================================================================
    % RMSE Frobenius: vs SNR (por ISR), vs ISR (por SNR)
    % =====================================================================
    for iISR_fix = 1:nISR
        fig_h = figure('Name', sprintf('RMSE Frob. vs SNR (ISR=%+d dB, init=%s)', ...
                       range_ISR_dB(iISR_fix), init_tag), ...
                       'NumberTitle','off', 'Position',[100 100 900 600]);
        hold on; grid on;
        for ic = 1:nCoupling
            coup_code = range_coupling(ic);
            if coup_code < 2, continue; end
            curve = rmse_eps_per_cen{ic}(:, iISR_fix);
            if all(isnan(curve)), continue; end
            curve_plot = max(curve, EPS_FLOOR);
            plot(range_SNR_dB, curve_plot, ...
                'LineStyle', styles_cen{coup_code+1}, ...
                'Marker', markers_cen{coup_code+1}, ...
                'Color', cmap_cen(coup_code+1, :), ...
                'LineWidth', 1.8, 'MarkerSize', 8, ...
                'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'DisplayName', short_labels_cen{coup_code+1});
        end
        set(gca, 'YScale', 'log');
        yline(EPS_FLOOR, ':', sprintf('piso visual = %.0e', EPS_FLOOR), ...
              'Color',[0.5 0.5 0.5], 'LabelHorizontalAlignment','left');
        xlabel('SNR (dB)'); ylabel('RMSE de ||C_{hat} - C_{true}||_F / ||C_{true}||_F');
        title(sprintf('RMSE Frobenius vs SNR  |  ISR = %+d dB  |  init: %s', ...
              range_ISR_dB(iISR_fix), init_tag));
        legend('Location','best','NumColumns',2);

        exportgraphics(fig_h, fullfile(outDir, ...
            sprintf('rmse_frobenius_vs_SNR_ISR_%+d_init_%s.png', range_ISR_dB(iISR_fix), init_tag)), ...
            'Resolution', 200);
    end

    if plot_vs_isr
    for iSNR_fix = 1:nSNR
        fig_h = figure('Name', sprintf('RMSE Frob. vs ISR (SNR=%+d dB, init=%s)', ...
                       range_SNR_dB(iSNR_fix), init_tag), ...
                       'NumberTitle','off', 'Position',[100 100 900 600]);
        hold on; grid on;
        for ic = 1:nCoupling
            coup_code = range_coupling(ic);
            if coup_code < 2, continue; end
            curve = rmse_eps_per_cen{ic}(iSNR_fix, :);
            if all(isnan(curve)), continue; end
            curve_plot = max(curve, EPS_FLOOR);
            plot(range_ISR_dB, curve_plot, ...
                'LineStyle', styles_cen{coup_code+1}, ...
                'Marker', markers_cen{coup_code+1}, ...
                'Color', cmap_cen(coup_code+1, :), ...
                'LineWidth', 1.8, 'MarkerSize', 8, ...
                'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'DisplayName', short_labels_cen{coup_code+1});
        end
        set(gca, 'YScale', 'log');
        yline(EPS_FLOOR, ':', sprintf('piso visual = %.0e', EPS_FLOOR), ...
              'Color',[0.5 0.5 0.5], 'LabelHorizontalAlignment','left');
        xlabel('ISR (dB)'); ylabel('RMSE de ||C_{hat} - C_{true}||_F / ||C_{true}||_F');
        title(sprintf('RMSE Frobenius vs ISR  |  SNR = %+d dB  |  init: %s', ...
              range_SNR_dB(iSNR_fix), init_tag));
        legend('Location','best','NumColumns',2);

        exportgraphics(fig_h, fullfile(outDir, ...
            sprintf('rmse_frobenius_vs_ISR_SNR_%+d_init_%s.png', range_SNR_dB(iSNR_fix), init_tag)), ...
            'Resolution', 200);
    end
    end   % --- if plot_vs_isr ---

    % =====================================================================
    % BER e EVM medios (apos beamforming): vs SNR (por ISR) e vs ISR (por SNR)
    % Duas curvas por cenario: uma para DAS, outra para MVDR/Capon.
    % =====================================================================
    % Pre-calcula matrizes (SNR x ISR) por cenario, agregando sobre angulos
    mean_BER_DAS_per_cen   = cell(nCoupling, 1);
    mean_BER_MVDR_per_cen  = cell(nCoupling, 1);
    mean_EVM_DAS_per_cen   = cell(nCoupling, 1);
    mean_EVM_MVDR_per_cen  = cell(nCoupling, 1);
    for ic_m = 1:nCoupling
        slab = agg_BER_DAS(ic_m, iInitPlot_agg, :, :, :, :);
        mean_BER_DAS_per_cen{ic_m}  = reshape(mean(slab, [5 6], 'omitnan'), [nSNR, nISR]);

        slab = agg_BER_MVDR(ic_m, iInitPlot_agg, :, :, :, :);
        mean_BER_MVDR_per_cen{ic_m} = reshape(mean(slab, [5 6], 'omitnan'), [nSNR, nISR]);

        slab = agg_EVM_DAS(ic_m, iInitPlot_agg, :, :, :, :);
        mean_EVM_DAS_per_cen{ic_m}  = reshape(mean(slab, [5 6], 'omitnan'), [nSNR, nISR]);

        slab = agg_EVM_MVDR(ic_m, iInitPlot_agg, :, :, :, :);
        mean_EVM_MVDR_per_cen{ic_m} = reshape(mean(slab, [5 6], 'omitnan'), [nSNR, nISR]);
    end

    % --- Funcao auxiliar inline para plotar uma metrica ---
    plot_metric = @(metric_per_cen, metric_name, file_prefix, ylabel_text, use_log) ...
        plot_metric_helper(metric_per_cen, metric_name, file_prefix, ylabel_text, ...
                           use_log, nCoupling, range_coupling, nISR, nSNR, ...
                           range_ISR_dB, range_SNR_dB, init_tag, ...
                           styles_cen, markers_cen, cmap_cen, short_labels_cen, ...
                           outDir, plot_vs_snr, plot_vs_isr);

    % BER do DAS
    plot_metric(mean_BER_DAS_per_cen,  'BER DAS',  'ber_das',  'BER medio (%)', false);
    % BER do MVDR/Capon
    plot_metric(mean_BER_MVDR_per_cen, 'BER MVDR', 'ber_mvdr', 'BER medio (%)', false);
    % EVM do DAS
    plot_metric(mean_EVM_DAS_per_cen,  'EVM DAS',  'evm_das',  'EVM medio (dB)', false);
    % EVM do MVDR/Capon
    plot_metric(mean_EVM_MVDR_per_cen, 'EVM MVDR', 'evm_mvdr', 'EVM medio (dB)', false);

    % =====================================================================
    % ACURACIA AGREGADA: um ponto por cenario (agregado sobre tudo)
    % =====================================================================
    rmse_phi_global = zeros(nCoupling, 1);
    rmse_eps_global = zeros(nCoupling, 1);
    for ic = 1:nCoupling
        rmse_phi_global(ic) = sqrt(mean(rmse_doa_per_cen{ic}(:).^2, 'omitnan'));
        rmse_eps_global(ic) = sqrt(mean(rmse_eps_per_cen{ic}(:).^2, 'omitnan'));
    end

    fig_h = figure('Name', sprintf('Acuracia agregada (init=%s)', init_tag), ...
                   'NumberTitle','off', 'Position',[100 100 900 700]);
    hold on; grid on;
    for ic = 1:nCoupling
        coup_code = range_coupling(ic);
        if coup_code < 2 || isnan(rmse_eps_global(ic)), continue; end
        method_label = doa_methods{method_idx(ic)};
        scatter(rmse_phi_global(ic), max(rmse_eps_global(ic), 1e-6), ...
                180, markers_cen{coup_code+1}, 'filled', ...
                'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
                'DisplayName', sprintf('%s (%s)', short_labels_cen{coup_code+1}, method_label));
    end
    xline(1, 'k:', '1° de erro de DoA');
    yline(0.1, 'k:', '\epsilon_F = 0.1');
    yline(1, 'r--', '\epsilon_F = 1 (sem comp.)');
    set(gca, 'YScale', 'log');
    xlabel('RMSE de \phi (graus)  [agregado]');
    ylabel('RMSE de ||C_{hat} - C_{true}||_F / ||C_{true}||_F  [agregado]');
    title(sprintf('Acuracia agregada: erro de DoA vs erro de C (init=%s)', init_tag));
    legend('Location','best');

    exportgraphics(fig_h, fullfile(outDir, ...
        sprintf('accuracy_agregada_init_%s.png', init_tag)), ...
        'Resolution', 200);
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

function plot_metric_helper(metric_per_cen, metric_name, file_prefix, ...
                            ylabel_text, use_log, ...
                            nCoupling, range_coupling, nISR, nSNR, ...
                            range_ISR_dB, range_SNR_dB, init_tag, ...
                            styles_cen, markers_cen, cmap_cen, ...
                            short_labels_cen, outDir, pv_snr, pv_isr)
%PLOT_METRIC_HELPER  Gera plots agregados de uma metrica (BER ou EVM)
%   vs SNR (uma figura por valor de ISR) e vs ISR (uma figura por SNR).
%   pv_snr/pv_isr controlam quais conjuntos sao gerados.

    % --- vs SNR (uma figura por ISR) ---
    if pv_snr
    for iISR_fix = 1:nISR
        fig_h = figure('Name', sprintf('%s vs SNR (ISR=%+d dB, init=%s)', ...
                       metric_name, range_ISR_dB(iISR_fix), init_tag), ...
                       'NumberTitle','off', 'Position',[100 100 900 600]);
        hold on; grid on;
        for ic = 1:nCoupling
            coup_code = range_coupling(ic);
            curve = metric_per_cen{ic}(:, iISR_fix);
            if all(isnan(curve)), continue; end
            plot(range_SNR_dB, curve, ...
                'LineStyle', styles_cen{coup_code+1}, ...
                'Marker', markers_cen{coup_code+1}, ...
                'Color', cmap_cen(coup_code+1, :), ...
                'LineWidth', 1.8, 'MarkerSize', 8, ...
                'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'DisplayName', short_labels_cen{coup_code+1});
        end
        if use_log, set(gca, 'YScale', 'log'); end
        xlabel('SNR (dB)'); ylabel(ylabel_text);
        title(sprintf('%s vs SNR  |  ISR = %+d dB  |  init: %s', ...
              metric_name, range_ISR_dB(iISR_fix), init_tag));
        legend('Location','best','NumColumns',2);

        exportgraphics(fig_h, fullfile(outDir, ...
            sprintf('%s_vs_SNR_ISR_%+d_init_%s.png', file_prefix, ...
                    range_ISR_dB(iISR_fix), init_tag)), ...
            'Resolution', 200);
    end
    end   % --- if pv_snr ---

    % --- vs ISR (uma figura por SNR) ---
    if pv_isr
    for iSNR_fix = 1:nSNR
        fig_h = figure('Name', sprintf('%s vs ISR (SNR=%+d dB, init=%s)', ...
                       metric_name, range_SNR_dB(iSNR_fix), init_tag), ...
                       'NumberTitle','off', 'Position',[100 100 900 600]);
        hold on; grid on;
        for ic = 1:nCoupling
            coup_code = range_coupling(ic);
            curve = metric_per_cen{ic}(iSNR_fix, :);
            if all(isnan(curve)), continue; end
            plot(range_ISR_dB, curve, ...
                'LineStyle', styles_cen{coup_code+1}, ...
                'Marker', markers_cen{coup_code+1}, ...
                'Color', cmap_cen(coup_code+1, :), ...
                'LineWidth', 1.8, 'MarkerSize', 8, ...
                'MarkerFaceColor', cmap_cen(coup_code+1, :), ...
                'DisplayName', short_labels_cen{coup_code+1});
        end
        if use_log, set(gca, 'YScale', 'log'); end
        xlabel('ISR (dB)'); ylabel(ylabel_text);
        title(sprintf('%s vs ISR  |  SNR = %+d dB  |  init: %s', ...
              metric_name, range_SNR_dB(iSNR_fix), init_tag));
        legend('Location','best','NumColumns',2);

        exportgraphics(fig_h, fullfile(outDir, ...
            sprintf('%s_vs_ISR_SNR_%+d_init_%s.png', file_prefix, ...
                    range_SNR_dB(iSNR_fix), init_tag)), ...
            'Resolution', 200);
    end
    end   % --- if pv_isr ---
end
