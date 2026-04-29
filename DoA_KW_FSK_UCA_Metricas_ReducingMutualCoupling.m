% Antes esse código era phase-only, agora ele é delay+phase
clear; clc;
%rng(1);
close all;
%% Parâmetros do array
M      = 8;           % nº de elementos do ULA
fc     = 2400e6; %500e6;         % Hz
c      = 3e8;
lambda = c/fc;
r = 0.25 * lambda;    % raio 1/4 λ

theta_sig_deg = 90;   % plano XY
theta_int_deg = 90;
sigma_erro_phi = 0; %10;

%nTrials  = 10; %360;   % nº de realizações
fs     = 288000;       % taxa de amostragem (Hz) para formar snapshots
N      = 21000;        % nº de amostras
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
range_ISR_dB = -20; %-6:1:-3;
range_snapshots = 2000:3000:N_DOA;
range_radius = [0.5]; %[0.25 0.2 0.15 0.1];
range_coupling = [1]; % com e sem acoplamento
range_phi = [75 125]; %-180:18:180;
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

C_true = eye(M);
for i = 1:M
    for j = 1:M
        d = mod(i-j, M);
        d = min(d, M-d);   % distância circular
        if d == 1
            c1_true = Z(i,j);
        elseif d == 2
            c2_true = Z(i,j);
        elseif i == j
            c3_true = Z(i,j);
        else
            c4_true = Z(i,j);
        end
    end
end


Coupling_matrices = zeros(M, M, nRadius); 
% Matrizes a ser usada nas simulações
for iRadius = 1:nRadius
    radius = range_radius(iRadius)*lambda;
    Ctx = Z;
    % Ctx = compute_Ctx_for_R(fc, M, radius, Z0);
    Coupling_matrices(:,:,iRadius) = Ctx; %inv(Ctx);
end

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

                            % Aplicando o Mutual Coupling
                            % Coupling_matrix = utils.coupling_matrix_uca(M, r, lambda, coupling_factor);
                            if Coupling == 0
                                Coupling_matrix = eye(M);
                                Coupling_matrix_sig = eye(M);
                                Coupling_matrix_int = eye(M);
                            else
                                % Ctx = compute_Ctx_for_R(fc, M, radius, Z0);
                                % Coupling_matrix = 1\Ctx;
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

                            if Coupling == 1
                                name_string = 'Coupling';
                            else
                                name_string = 'No Coupling';
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

                            % ----- ESTIMAÇÃO DA MCM -----
                            if abs(phi_sig_deg - 75) < 1e-9
                                a_sig = utils.steering_vec_uca(M, radius, lambda, theta_sig_deg, phi_sig_deg);  % vetor de direção Mx1

                                b_hat = X * q / (q' * q);   % 4 x 1

                                % Etapa 2: ajustar c1 e c2

                                % chute inicial
                                p0 = zeros(1,8);

                                cost_fun = @(p) cost_mcm_circulant4(p, a, b_hat);

                                opts = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10,...
                                    'MaxIter',5000,'MaxFunEvals',20000);

                                p_est = fminsearch(cost_fun, p0, opts);

                                c1_est = p_est(1) + 1j*p_est(2);
                                c2_est = p_est(3) + 1j*p_est(4);
                                c3_est = p_est(5) + 1j*p_est(6);
                                c4_est = p_est(7) + 1j*p_est(8);

                                C_est = build_C_circulant_uca4(c1_est, c2_est, c3_est, c4_est);

                                % Resultados
                                disp('=== Coeficientes verdadeiros ===')
                                disp(['c1_true = ', num2str(c1_true)])
                                disp(['c2_true = ', num2str(c2_true)])
                                disp(['c3_true = ', num2str(c3_true)])
                                disp(['c4_true = ', num2str(c4_true)])

                                disp('=== Coeficientes estimados ===')
                                disp(['c1_est  = ', num2str(c1_est)])
                                disp(['c2_est  = ', num2str(c2_est)])
                                disp(['c3_est  = ', num2str(c3_est)])
                                disp(['c4_est  = ', num2str(c4_est)])

                                fprintf('Erro relativo Frobenius = %.6e\n', ...
                                    norm(C_est - C_true, 'fro') / norm(C_true, 'fro'));

                                % Comparações
                                b_model_true = Coupling_matrix * a_sig;
                                b_model_est  = C_est  * a_sig;

                                % remove ambiguidade escalar para comparar vetores
                                alpha_true = (b_model_true' * b_hat) / (b_model_true' * b_model_true);
                                alpha_est  = (b_model_est'  * b_hat) / (b_model_est'  * b_model_est);

                                figure;
                                subplot(1,2,1);
                                imagesc(abs(C_true)); colorbar; axis equal tight;
                                title('|C true|');

                                subplot(1,2,2);
                                imagesc(abs(C_est)); colorbar; axis equal tight;
                                title('|C est|');

                                figure;
                                subplot(2,1,1);
                                plot(1:M, abs(b_hat), 'ko-','LineWidth',1.5); hold on;
                                plot(1:M, abs(alpha_true*b_model_true), 'bx--','LineWidth',1.5);
                                plot(1:M, abs(alpha_est*b_model_est), 'r*-','LineWidth',1.5);
                                grid on;
                                xlabel('Indice da antena');
                                ylabel('Magnitude');
                                legend('|b_{hat}|','|b_{true}| ajustado','|b_{est}| ajustado','Location','best');
                                title('Comparacao do vetor espacial efetivo - magnitude');

                                subplot(2,1,2);
                                plot(1:M, unwrap(angle(b_hat)), 'ko-','LineWidth',1.5); hold on;
                                plot(1:M, unwrap(angle(alpha_true*b_model_true)), 'bx--','LineWidth',1.5);
                                plot(1:M, unwrap(angle(alpha_est*b_model_est)), 'r*-','LineWidth',1.5);
                                grid on;
                                xlabel('Indice da antena');
                                ylabel('Fase (rad)');
                                legend('angle(b_{hat})','angle(b_{true}) ajustado','angle(b_{est}) ajustado','Location','best');
                                title('Comparacao do vetor espacial efetivo - fase');
                            end

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
                                    figure('Position', [50 50 1500 900], 'Color', 'w', ...
                                        'Name', fig_name);
            
                                    % (1) Sinal ideal (sem acoplamento)
                                    subplot(5,1,1);
                                    plot(real(q(1:N_plot)), 'b-', 'LineWidth', 0.8); hold on;
                                    plot(imag(q(1:N_plot)), 'r-', 'LineWidth', 0.8);
                                    grid on;
                                    ylabel('Amplitude');
                                    title(sprintf('(a) Sem ideal', ant_idx));
                                    % legend('Re\{x_{ideal}(t)\}', 'Simbolo', 'Location', 'northeast');

                                    % (1) Sinal ideal (sem acoplamento)
                                    subplot(5,1,2);
                                    plot(real(Xq(ant_idx,1:N_plot)), 'b-', 'LineWidth', 0.8); hold on;
                                    plot(imag(Xq(ant_idx,1:N_plot)), 'r-', 'LineWidth', 0.8);
                                    grid on;
                                    ylabel('Amplitude');
                                    title(sprintf('(b) Sem acoplamento - Antena %d', ant_idx));
                                    % legend('Re\{x_{ideal}(t)\}', 'Simbolo', 'Location', 'northeast');

                                    % (2) Sinal com acoplamento (sem beamforming)
                                    subplot(5,1,3);
                                    plot(real(X(ant_idx,1:N_plot)), 'b-', 'LineWidth', 0.8); hold on;
                                    plot(imag(X(ant_idx,1:N_plot)), 'r-', 'LineWidth', 0.8);
                                    grid on;
                                    ylabel('Amplitude');
                                    title(sprintf('(b) Com acoplamento (C^{-1}) - Antena %d', ant_idx));
                                    % legend('Re\{x_{acoplado}(t)\}', 'Simbolo', 'Location', 'northeast');

                                    % (3) Saida Capon (MVDR)
                                    subplot(5,1,4);
                                    plot(real(y_capon(1:N_plot)), 'b-', 'LineWidth', 0.8); hold on;
                                    plot(imag(y_capon(1:N_plot)), 'r-', 'LineWidth', 0.8);
                                    grid on;
                                    ylabel('Amplitude');
                                    title('(c) Saida Capon (MVDR) - apos acoplamento + beamforming');
                                    % legend('Re\{y_{Capon}(t)\}', 'Simbolo', 'Location', 'northeast');

                                    % (4) Saida DAS
                                    subplot(5,1,5);
                                    plot(real(y_das(1:N_plot)), 'b-', 'LineWidth', 0.8); hold on;
                                    plot(imag(y_das(1:N_plot)), 'r-', 'LineWidth', 0.8);
                                    grid on;
                                    hold off; grid on;
                                    xlabel('Tempo (ms)');
                                    ylabel('Amplitude');
                                    title('(d) Saida DAS - apos acoplamento + beamforming');
                                    % legend('Re\{y_{DAS}(t)\}', 'Simbolo', 'Location', 'northeast');

                                    sgtitle(sprintf('Cadeia de sinal | Coupling %d | SNR=%d dB | ISR=%d dB | R=%.2f\\lambda | \\phi_{sig}=%.0f° | \\phi_{int}=%.0f°', ...
                                        Coupling, SNR_dB, ISR_dB, range_radius(iRadius), phi_sig_deg, phi_int_deg), ...
                                        'FontSize', 13, 'FontWeight', 'bold');

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

                                if Coupling == 1
                                    name_string = 'Coupling';
                                else
                                    name_string = 'No Coupling';
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