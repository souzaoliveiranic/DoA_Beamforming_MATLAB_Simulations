% Antes esse código era phase-only, agora ele é delay+phase
clear; clc; 
%rng(1);
close all;
%% Parâmetros do array
M      = 8;           % nº de elementos do ULA
fc     = 500e6;         % Hz 
c      = 3e8;        
lambda = c/fc;        
r = 0.25 * lambda;    % raio 1/4 λ
% coupling_factor = 0.3;

sigma_erro_phi = 0; %10;

%nTrials  = 10; %360;   % nº de realizações
fs     = 288000;       % taxa de amostragem (Hz) para formar snapshots
N      = 21000;        % nº de amostras
N_DOA  = 3100;
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
% Modo de varredura angular (2D - azimute e elevacao):
%   'cross'    -> phi e theta varrem grades via produto cartesiano.
%                 nPhi*nTheta * nPhi*nTheta simulacoes por ponto operacional.
%   'aleatory' -> (phi_sig, theta_sig) e (phi_int, theta_int) sorteados
%                 aleatoriamente, respeitando separacao angular esferica
%                 minima. n_rand_angles pares por ponto operacional.
% =========================================================================
experiment_mode = 'aleatory';   % 'cross' | 'aleatory'

% Parametros do modo 'aleatory'
n_rand_angles  = 500;            % quantos pares sortear por ponto operacional
min_sep_deg    = 10;            % separacao angular esferica minima (graus)
theta_min_rand = 10;            % limite inferior de theta no sorteio (evita zenith)
theta_max_rand = 90;            % limite superior de theta no sorteio

range_SNR_dB    = -6:3:6;
range_ISR_dB    = -12:3:6;
range_snapshots = N_DOA; %2000:3000:N_DOA;
range_radius    = [0.22];
range_coupling  = [0 1];        % 0 = sem coupling, 1 = com coupling

% Grades angulares verdadeiras (usadas em 'cross')
range_phi_true   = -180:18:180;       % azimutes verdadeiros
range_theta_true = 0:15:90;           % elevacoes verdadeiras [0,90] step 15

switch experiment_mode
    case 'cross'
        % Produto cartesiano de (phi, theta) para sig e int
        [PhiG, ThetaG] = meshgrid(range_phi_true, range_theta_true);
        sig_dirs = [PhiG(:), ThetaG(:)];  % nPairs x 2 (col1=phi, col2=theta)
        int_dirs = sig_dirs;              % mesmo conjunto para interferente
        % Pares serao formados via duplo loop excluindo iguais
        nSigDirs = size(sig_dirs, 1);
        nIntDirs = size(int_dirs, 1);
        nPairsPerOp = nSigDirs * nIntDirs;   % cap superior; pares iguais sao puladas
        pair_mode = 'cross';

    case 'aleatory'
        % Sorteia n_rand_angles pares (sig, int) respeitando separacao
        % esferica minima. Angulos FIXOS para todos os pontos operacionais
        % (so SNR/ISR/ruido variam entre pontos).
        rng(12345, 'twister');
        sig_dirs = zeros(n_rand_angles, 2);   % [phi, theta]
        int_dirs = zeros(n_rand_angles, 2);
        for kk = 1:n_rand_angles
            % Sorteia SOI
            ps = -180 + 360*rand;
            ts = theta_min_rand + (theta_max_rand - theta_min_rand)*rand;
            % Sorteia interferente com separacao esferica minima
            while true
                pi_ = -180 + 360*rand;
                ti_ = theta_min_rand + (theta_max_rand - theta_min_rand)*rand;
                d_sph = spherical_angular_distance(ts, ps, ti_, pi_);
                if d_sph >= min_sep_deg
                    break;
                end
            end
            sig_dirs(kk, :) = [ps, ts];
            int_dirs(kk, :) = [pi_, ti_];
        end
        nPairsPerOp = n_rand_angles;
        pair_mode = 'paired';

    otherwise
        error('experiment_mode invalido: %s', experiment_mode);
end

nMethods = 4;
methods = ["KW","DAS","MPDR","MUSIC"];

nSNR        = numel(range_SNR_dB);
nISR        = numel(range_ISR_dB);
nRadius     = numel(range_radius);
nSnapshots  = numel(range_snapshots);
nCoupling   = numel(range_coupling);

% =========================================================================
% Grade de busca 2D para metodos classicos (DAS, MPDR, MUSIC)
% theta passo 1 grau, phi passo 0.5 grau (sweep refinado)
% =========================================================================
theta_scan_step = 1;
phi_scan_step   = 0.5;
theta_scan_grid = 1:theta_scan_step:90;          % evita zenith (degenerescencia)
phi_scan_grid   = -180:phi_scan_step:180;
nThetaScan = numel(theta_scan_grid);
nPhiScan   = numel(phi_scan_grid);

% Vetores planos (theta, phi) para cada ponto da grade 2D
[PhiScanG, ThetaScanG] = meshgrid(phi_scan_grid, theta_scan_grid);
theta_scan_flat = ThetaScanG(:);   % nGridPts x 1
phi_scan_flat   = PhiScanG(:);     % nGridPts x 1
nGridPts        = numel(theta_scan_flat);
fprintf('Grade de busca 2D: %d pontos (theta=%d, phi=%d)\n', ...
        nGridPts, nThetaScan, nPhiScan);

TotalSim = nSNR*nISR*nRadius*nCoupling*nSnapshots*nPairsPerOp;
iTotal   = 0;

% Armazenamento (ultima dim agora indexa o par angular dentro do ponto operacional)
RMSE = zeros(nMethods, nSNR, nISR, nRadius, nCoupling, nSnapshots, nPairsPerOp);

%% ======= Coupling Matrix Estimation =======

Z0 = 50; % Impedância de referência

% R_list = [0.05 0.10 0.15 0.20 0.30];   % metros (ajuste)
R_list = [0.01 0.025 0.125 0.25 0.375 0.5];   % metros (ajuste)
use_dB = false;

% --------- Pré-alocação ---------
Ccols = cell(numel(R_list),1);   % guarda |C(:,1)| para cada R
Dists = cell(numel(R_list),1);   % guarda distâncias até antena 1

% --------- Loop nos raios ---------
for k = 1:numel(R_list)
    R = R_list(k)*lambda;

    % 1) Calcula Ctx para esse raio (via S->Z e fórmula da pag 29)
    Ctx = compute_Ctx_for_R(fc, M, R, Z0);
    Ctx = 1\Ctx;
    % 2) Coluna 1
    c1 = Ctx(:,1);

    if use_dB
        cplot = 20*log10(abs(c1)+eps);
    else
        cplot = abs(c1);
    end
    Ccols{k} = cplot;

    % 3) Distâncias do UCA até a antena 1
    phi = phi0 + (0:M-1).' * (2*pi/M);
    P = [R*cos(phi), R*sin(phi), zeros(M,1)];  % Mx3

    d = zeros(M,1);
    for i = 1:M
        d(i) = norm(P(i,:) - P(1,:));
    end
    Dists{k} = d;
end

% --------- Plot 1: |C(i,1)| vs índice ---------
figure; hold on;
x = 1:M;
for k = 1:numel(R_list)
    plot(x, Ccols{k}, 'o-','LineWidth',1.3);
end
grid on;
xlabel('Índice da antena i');
if use_dB
    ylabel("|C_{i,1}| (dB)");
else
    ylabel("|C_{i,1}|");
end
 
title('Comparação da Magnitude do Acoplamento para diferentes raios (r) do UCA');
% legend(compose('r = %.3f m', R_list), 'Location','best');
legend(compose('R = %.3f $\\lambda_0$', R_list), ...
       'Interpreter','latex', ...
       'Location','best');

% --------- Plot 2: |C(i,1)| vs distância até antena 1 ---------
figure; hold on;
for k = 1:numel(R_list)
    d = Dists{k};
    cplot = Ccols{k};

    % ordenar por distância para o gráfico ficar "bonito"
    [ds, ord] = sort(d);
    cs = cplot(ord);

    plot(ds, cs, 'o-','LineWidth',1.3);
end
grid on;
xlabel('Distância até a antena 1 (m)');
if use_dB
    ylabel("|C_{i,1}| (dB)");
else
    ylabel("|C_{i,1}|");
end
title('Acoplamento vs distância (referência: antena 1) para diferentes R');
% legend(compose('R = %.3f \lambda', R_list), 'Location','best');
legend(compose('R = %.3f $\\lambda_0$', R_list), ...
       'Interpreter','latex', ...
       'Location','best');

Coupling_matrices = zeros(M, M, nRadius); 
% Matrizes a ser usada nas simulações
for iRadius = 1:nRadius
    radius = range_radius(iRadius)*lambda;
    Ctx = compute_Ctx_for_R(fc, M, radius, Z0);
    Coupling_matrices(:,:,iRadius) = Ctx;%inv(Ctx);
end

%% ======= Pré-computação dos steering vectors da grade 2D =======
% A matriz de steering vectors da grade de busca (theta, phi) e' fixa por raio.
% Pre-computar fora do loop evita milhoes de chamadas a steering_vec_uca.
fprintf('Pre-computando matriz de steering vectors 2D (%d pontos x %d raios)...\n', ...
        nGridPts, nRadius);
tic_pre = tic;
A_scan_cell = cell(nRadius, 1);
for iRadius_pre = 1:nRadius
    radius_pre = range_radius(iRadius_pre) * lambda;
    A = zeros(M, nGridPts);
    for kk = 1:nGridPts
        A(:, kk) = utils.steering_vec_uca(M, radius_pre, lambda, ...
                                          theta_scan_flat(kk), phi_scan_flat(kk));
    end
    A_scan_cell{iRadius_pre} = A;
end
fprintf('  ...pre-computacao concluida em %.2f s\n', toc(tic_pre));

% Vetor unitario das direcoes verdadeiras para metrica esferica
% (sera (re)usado dentro do loop)

%% ======= DoA KW vs Delay and Sum vs Capon =======

sweep_t0 = tic;
fprintf('\n=== INICIO DO LOOP PRINCIPAL ===\n');
fprintf('Total de simulacoes: %d  (modo %s, %d pares angulares)\n', ...
        TotalSim, experiment_mode, nPairsPerOp);
fprintf('================================\n\n');

for iSNR = 1:nSNR
    for iISR = 1:nISR
        for iRadius = 1:nRadius
            for iCoupling = 1:nCoupling

                SNR_dB = range_SNR_dB(iSNR);
                ISR_dB = range_ISR_dB(iISR);
                radius = range_radius(iRadius)*lambda;
                Coupling = range_coupling(iCoupling);

                % Steering scan matrix (pre-computada) para este raio
                A_scan = A_scan_cell{iRadius};

                % =========================================================
                % Geracao da lista de pares angulares para este ponto operacional
                % =========================================================
                if strcmp(pair_mode, 'cross')
                    % Produto cartesiano (sig_dirs x int_dirs), descartando iguais
                    pair_list = zeros(nPairsPerOp, 4);  % [phi_s, theta_s, phi_i, theta_i]
                    cnt = 0;
                    for is = 1:size(sig_dirs, 1)
                        for ii_ = 1:size(int_dirs, 1)
                            if is == ii_  % mesma direcao -> pula
                                continue;
                            end
                            cnt = cnt + 1;
                            pair_list(cnt, :) = [sig_dirs(is, 1), sig_dirs(is, 2), ...
                                                 int_dirs(ii_, 1), int_dirs(ii_, 2)];
                        end
                    end
                    pair_list = pair_list(1:cnt, :);
                else  % 'paired' (modo 'aleatory')
                    pair_list = [sig_dirs, int_dirs];   % colunas [phi_s,theta_s,phi_i,theta_i]
                end
                nPairs_actual = size(pair_list, 1);

                ii = 1;   % indice linear do par para armazenamento
                for iPair = 1:nPairs_actual
                    phi_sig_deg   = pair_list(iPair, 1);
                    theta_sig_deg = pair_list(iPair, 2);
                    phi_int_deg   = pair_list(iPair, 3);
                    theta_int_deg = pair_list(iPair, 4);

                    j = 1;

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
                        X = Coupling_matrix * Xsig + Coupling_matrix * Xint + Xn;

                        for K = range_snapshots
                            iTotal = iTotal + 1;
                            % --- Progresso da simulacao ---
                            elapsed = toc(sweep_t0);
                            if iTotal > 1
                                eta_sec = elapsed * (TotalSim - iTotal) / (iTotal - 1);
                                eta_str = datestr(seconds(eta_sec), 'HH:MM:SS');
                            else
                                eta_str = '--:--:--';
                            end
                            fprintf(['--> Sim %d / %d (%.1f%%)  |  SNR=%+d ISR=%+d Coup=%d r=%d  |  ' ...
                                     'sig=(%+6.1f, %4.1f)  int=(%+6.1f, %4.1f)  K=%d  |  ' ...
                                     'decorrido %s  ETA %s\n'], ...
                                    iTotal, TotalSim, 100*iTotal/TotalSim, ...
                                    SNR_dB, ISR_dB, Coupling, iRadius, ...
                                    phi_sig_deg, theta_sig_deg, phi_int_deg, theta_int_deg, K, ...
                                    datestr(seconds(elapsed), 'HH:MM:SS'), eta_str);

                            % ----- DoA KW (2D nativo) -----
                            beta = 2*pi*(0:M-1)'/M;
                            [theta_hat_KW, phi_hat_KW] = doa_kw_uca(X(:,1:K), q(1:K).', radius, lambda, beta);

                            % Matriz de covariância
                            Rxx = (X(:,1:K)*X(:,1:K)')/size(X(:,1:K),2);
                            delta = 1e-3 * trace(Rxx)/M;
                            Rxx_dl = Rxx + delta*eye(M);

                            % Decomposição espectral
                            [eigvec, eigval] = eig(Rxx_dl);
                            [~, idx_eig] = sort(diag(eigval), 'descend');
                            E = eigvec(:, idx_eig);

                            % Número de fontes conhecidas
                            Ksrc = 2;
                            En = E(:, Ksrc+1:end);
                            EnEnH = En * En';

                            % ===== Varredura 2D (theta, phi) VETORIZADA =====
                            % Para A_scan (M x nGridPts) e B (M x M) hermitiana,
                            % o vetor [a_k' * B * a_k] (k=1..nGridPts) e' sum(conj(A) .* (B*A), 1).
                            Rinv = inv(Rxx_dl);
                            BA_DAS   = Rxx   * A_scan;
                            BA_MVDR  = Rinv  * A_scan;
                            BA_MUSIC = EnEnH * A_scan;

                            P_DAS_flat   = abs(  sum(conj(A_scan) .* BA_DAS,   1) );
                            denom_mvdr   = real( sum(conj(A_scan) .* BA_MVDR,  1) );
                            P_MVDR_flat  = 1 ./ max(denom_mvdr, eps);
                            denom_music  = real( sum(conj(A_scan) .* BA_MUSIC, 1) );
                            P_MUSIC_flat = 1 ./ max(denom_music, eps);

                            % Indices dos picos (na grade 2D plana)
                            [~, i_das  ] = max(P_DAS_flat);
                            [~, i_mvdr ] = max(P_MVDR_flat);
                            [~, i_music] = max(P_MUSIC_flat);

                            % Recupera (theta, phi) estimados de cada metodo
                            theta_DAS   = theta_scan_flat(i_das);    phi_DAS   = phi_scan_flat(i_das);
                            theta_MVDR  = theta_scan_flat(i_mvdr);   phi_MVDR  = phi_scan_flat(i_mvdr);
                            theta_MUSIC = theta_scan_flat(i_music);  phi_MUSIC = phi_scan_flat(i_music);

                            % ===== Erro angular ESFERICO (geodesico) =====
                            err_KW    = spherical_angular_distance(theta_sig_deg, phi_sig_deg, theta_hat_KW,  phi_hat_KW);
                            err_DAS   = spherical_angular_distance(theta_sig_deg, phi_sig_deg, theta_DAS,     phi_DAS);
                            err_MVDR  = spherical_angular_distance(theta_sig_deg, phi_sig_deg, theta_MVDR,    phi_MVDR);
                            err_MUSIC = spherical_angular_distance(theta_sig_deg, phi_sig_deg, theta_MUSIC,   phi_MUSIC);

                            RMSE(:, iSNR, iISR, iRadius, iCoupling, j, ii) = ...
                                [err_KW; err_DAS; err_MVDR; err_MUSIC];
                            j = j + 1;
                        end
                        ii = ii + 1;
                end  % iPair
            end  % iCoupling
        end  % iRadius
    end  % iISR
end  % iSNR

fprintf('\nLoop principal concluido em %s\n', ...
        datestr(seconds(toc(sweep_t0)), 'HH:MM:SS'));

RMSE_mean = mean(RMSE, 7);   % média ao longo dos pares angulares
%K_fixed = range_snapshots(end);
K_fixed = numel(range_snapshots);




%% ======= Gráficos de SNR, ISR, CF e Snapshots =======

methods = ["KW","DAS","Capon","MUSIC"];
markers = ["o-","x-","s-","d-","^-","v-","*-","+-"];
colors =  ["red", "green", "blue", "black", "magenta", "cyan", "yellow"];
line_style =   ["-", "--", ":","-."]; % linha continua sem acoplamento % linha tracejada com acoplamento

% Criar pasta 'graficos'
outDir = fullfile(pwd, 'Novos Graficos EUSIPCO - SBRT');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%  GRÁFICOS: Varredura de SNR

for iRadius = 1:nRadius
    for iISR = 1:nISR

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
            ylim([-10 100]);
            % title(sprintf('RMSE vs SNR | ISR=%d dB | Snapshots=%d', ...
            %     range_ISR_dB(iISR), range_snapshots(end)));

            % legend(legend_strings, 'Location', 'northwest','Interpreter','latex');
            legend(legend_strings, 'Location', 'best','Interpreter','latex');

            fileName = sprintf('RMSE_vs_SNR_ISR_%ddB_Snap_%d.png', ...
                range_ISR_dB(iISR), range_snapshots(end));

            exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);

            fileName = sprintf('RMSE_vs_SNR_ISR_%ddB_Snap_%d.tex', ...
                range_ISR_dB(iISR), range_snapshots(end));

            %cleanfigure;
            % Define LaTeX macros for width and height (you will define these in your .tex file)
            matlab2tikz(fullfile(outDir, fileName), 'width', '\figurewidth', 'height', '\figureheight');
    end
end

%%  GRÁFICOS: Varredura de ISR

for iRadius = 1:nRadius
    for iSNR = 1:nSNR
        fig = figure('Name', "RMSE vs ISR | SNR=" + string(range_SNR_dB(iSNR)) + " (r = " + string(range_radius(iRadius)) + ")", ...
            'NumberTitle', 'off');
        hold on; grid on;

        legend_strings = strings(0);   % inicializa como array de strings vazio

        for iCoupling = 1:nCoupling
            for m = 1:nMethods
                y = squeeze(RMSE_mean(m, iSNR, :, iRadius, iCoupling, K_fixed));
                plot(range_ISR_dB, y, markers(m), 'Color', colors(m), 'LineWidth', 2, 'LineStyle', line_style(iCoupling));
            end
        end
        for m = 1:nMethods
            legend_strings(end+1) = methods(m);
        end

        xlabel("ISR (dB)",'Interpreter','latex');
        ylabel("RMSE (degrees)",'Interpreter','latex');
        yscale("log");
        ytickformat('%.2f') % 'f' forces fixed-point notation instead of scientific
        ylim([-10 100]);
        % title(sprintf('RMSE vs ISR | SNR=%d dB | Snapshots=%d', ...
        %     range_SNR_dB(iSNR), range_snapshots(end)));

        % legend(legend_strings, 'Location', 'northwest','Interpreter','latex');
        legend(legend_strings, 'Location', 'best','Interpreter','latex');

        fileName = sprintf('RMSE_vs_ISR_SNR_%ddB_Snap_%d.png', ...
            range_SNR_dB(iSNR), range_snapshots(end));

        exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);

        fileName = sprintf('RMSE_vs_ISR_SNR_%ddB_Snap_%d.tex', ...
            range_SNR_dB(iSNR), range_snapshots(end));

        %    cleanfigure;
        % Define LaTeX macros for width and height (you will define these in your .tex file)
        matlab2tikz(fullfile(outDir, fileName), 'width', '\figurewidth', 'height', '\figureheight');
    end
end

return;

%%  GRÁFICOS: Varredura de SNR por Raio

for m = 1:nMethods
    for iISR = 1:nISR
        if (range_ISR_dB(iISR) == -6)
            fig = figure('Name', "RMSE vs SNR | ISR=" + string(range_ISR_dB(iISR)) + " (" + methods(m) + ")", ...
                'NumberTitle', 'off');
            hold on; grid on;

            legend_strings = strings(0);   % inicializa como array de strings vazio

            for iCoupling = 1:nCoupling
                for iRadius = 3:nRadius
                    plot(range_SNR_dB, squeeze(RMSE_mean(m, :, iISR, iRadius, iCoupling, K_fixed)), ...
                        markers(iRadius), 'Color', colors(iRadius),'LineWidth', 2, 'LineStyle', line_style(iCoupling));
                end
            end
            for iRadius = 3:nRadius
                legend_strings(end+1) = "r = " + string(range_radius(iRadius));
            end

            xlabel("SNR (dB)",'Interpreter','latex');
            ylabel("RMSE (degrees)",'Interpreter','latex');
            yscale("log");
            ytickformat('%.2f') % 'f' forces fixed-point notation instead of scientific
            ylim([-10 100]);
            % title(sprintf('RMSE vs SNR | ISR=%d dB | Snapshots=%d', ...
            %     range_ISR_dB(iISR), range_snapshots(end)));

            % legend(legend_strings, 'Location', 'northwest','Interpreter','latex');
            legend(legend_strings, 'Location', 'best','Interpreter','latex');

            fileName = sprintf('RMSE_vs_SNR_ISR_%ddB_Snap_%d.png', ...
                range_ISR_dB(iISR), range_snapshots(end));

            exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);

            fileName = sprintf('RMSE_vs_SNR_ISR_%ddB_Snap_%d.tex', ...
                range_ISR_dB(iISR), range_snapshots(end));

            %cleanfigure;
            % Define LaTeX macros for width and height (you will define these in your .tex file)
            matlab2tikz(fullfile(outDir, fileName), 'width', '\figurewidth', 'height', '\figureheight');
        end
    end
end


for iISR = 1:nISR
    if (range_ISR_dB(iISR) == -6)
        fig = figure('Name', "RMSE vs SNR | ISR=" + string(range_ISR_dB(iISR)), ...
            'NumberTitle', 'off');
        hold on; grid on;

        legend_strings = strings(0);   % inicializa como array de strings vazio

        for m = 1:nMethods
            for iCoupling = 1:nCoupling
                if(range_coupling(iCoupling) == 1)
                    for iRadius = 3:nRadius
                        plot(range_SNR_dB, squeeze(RMSE_mean(m, :, iISR, iRadius, iCoupling, K_fixed)), ...
                            markers(m), 'Color', colors(iRadius),'LineWidth', 2, 'LineStyle', line_style(m));
                    end
                end
            end
        end
        for iRadius = 3:nRadius
            legend_strings(end+1) = "r = " + string(range_radius(iRadius));
        end

        xlabel("SNR (dB)",'Interpreter','latex');
        ylabel("RMSE (degrees)",'Interpreter','latex');
        yscale("log");
        ytickformat('%.2f') % 'f' forces fixed-point notation instead of scientific
        ylim([-10 100]);
        % title(sprintf('RMSE vs SNR | ISR=%d dB | Snapshots=%d', ...
        %     range_ISR_dB(iISR), range_snapshots(end)));

        % legend(legend_strings, 'Location', 'northwest','Interpreter','latex');
        legend(legend_strings, 'Location', 'best','Interpreter','latex');

        fileName = sprintf('RMSE_vs_SNR_ISR_%ddB_Snap_%d.png', ...
            range_ISR_dB(iISR), range_snapshots(end));

        exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);

        fileName = sprintf('RMSE_vs_SNR_ISR_%ddB_Snap_%d.tex', ...
            range_ISR_dB(iISR), range_snapshots(end));

        %cleanfigure;
        % Define LaTeX macros for width and height (you will define these in your .tex file)
        matlab2tikz(fullfile(outDir, fileName), 'width', '\figurewidth', 'height', '\figureheight');
    end
end




return;

%%  GRÁFICOS: Snapshots

% for iSNR = 1:nSNR
%     for iISR = 1:nISR
%         fig = figure('Name', "RMSE vs Snapshots | SNR=" + string(range_SNR_dB(iSNR)) + ...
%             " | ISR=" + string(range_ISR_dB(iISR)), ...
%             'NumberTitle', 'off');
%         hold on; grid on;
% 
%         legend_strings = strings(0);   % inicializa como array de strings vazio
% 
%         for iRadius = 1:nRadius
%             for iCoupling = 1:nCoupling
%                 for m = 1:nMethods
%                     y = squeeze(RMSE_mean(m, iSNR, iISR, iRadius, iCoupling, :));
%                     plot(range_snapshots, y, markers(iRadius), 'Color', colors(m), 'LineWidth', 2, 'LineStyle', line_style(iRadius));
%                 end
%             end
%         end
%         for m = 1:nMethods
%             legend_strings(end+1) = methods(m);
%         end
% 
%         xlabel("Snapshots (K)",'Interpreter','latex');
%         ylabel("RMSE (degrees)",'Interpreter','latex');
%         yscale("log");
%         ytickformat('%.2f') % 'f' forces fixed-point notation instead of scientific
%         % title(sprintf('RMSE vs Snapshots | SNR=%d dB | ISR=%d dB', ...
%         %     range_SNR_dB(iSNR), range_ISR_dB(iISR)));
% 
%         legend(legend_strings, 'Location', 'northeast','Interpreter','latex');
%         ylim([0 90]);
%         fileName = sprintf('RMSE_vs_Snap_SNR_%ddB_ISR_%ddB.png', ...
%             range_SNR_dB(iSNR), range_ISR_dB(iISR));
% 
%         exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);
% 
%         fileName = sprintf('RMSE_vs_Snap_SNR_%ddB_ISR_%ddB.tex', ...
%             range_SNR_dB(iSNR), range_ISR_dB(iISR));
% 
%         cleanfigure;
%         % Define LaTeX macros for width and height (you will define these in your .tex file)
%         matlab2tikz(fullfile(outDir, fileName), 'width', '\figurewidth', 'height', '\figureheight');
%     end
% end




%% ======= Calculo do Beamforming  =======

% Testando outros ranges no beamforming
range_SNR_dB = -9:3:3; %-9:3:9;
range_ISR_dB = -6:3:3; %-9:3:3;

nulo_no_interferidor = 0;

TotalSim = nSNR*nISR*nRadius*nCoupling*nPhi*nPhi;
iTotal = 0;
% Armazenamento
BER = zeros(4, nSNR, nISR, nRadius, nPhi*nPhi);
EVM = zeros(4, nSNR, nISR, nRadius, nPhi*nPhi);

for iSNR = 1:nSNR
    for iISR = 1:nISR
        for iRadius = 1:nRadius
            for iCoupling = 1:nCoupling

                SNR_dB = range_SNR_dB(iSNR);
                ISR_dB = range_ISR_dB(iISR);
                radius = range_radius(iRadius)*lambda;
                Coupling = range_coupling(iCoupling);

                ii = 1;
                for iPhi = 1:nPhi
                    phi_sig_deg = range_phi(iPhi);

                    for iiPhi = 1:nPhi
                        phi_int_deg = range_phi(iiPhi);

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
                            g_emb_sig = zeros(1,M);
                            g_emb_int = zeros(1,M);
                            Ctx = compute_Ctx_for_R(fc, M, radius, Z0);
                            Coupling_matrix = 1\Ctx;
                        end
                        X = Coupling_matrix * X;
                        % X = Coupling_matrix_sig * Xsig + Coupling_matrix_int * Xint + Xn;
                        % X = Coupling_matrix * Xsig + Coupling_matrix * Xint + Xn;


                        iTotal = iTotal + 1;
                        % ----- DoA KW
                        fprintf('--> Sim %d / %d , ii = %d , j = %d , K = %d / %d \n', iTotal, TotalSim, ...
                            ii, j, K, N_DOA);

                        % ----- BEAMFORMING -----

                        % ----- Restrições angulares

                        a_sig = utils.steering_vec_uca(M, radius, lambda, theta_sig_deg, phi_sig_deg);  % vetor de direção Mx1

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

                        % ======= Demodulação =======

                        [bits_hat_q, BER_q, pam_rx_mf0, sym_rx0]      = utils.fsk2_demod(qn, bits, Rs, sps, alpha, span, fd);
                        [bits_hat_controle, BER_controle, pam_rx_mf4, sym_rx4]      = utils.fsk2_demod(Xsig(1,:), bits, Rs, sps, alpha, span, fd);
                        [bits_hat_in, BER_in, pam_rx_mf1, sym_rx1]    = utils.fsk2_demod(X(1,:), bits, Rs, sps, alpha, span, fd);
                        [bits_hat_das, BER_das, pam_rx_mf2, sym_rx2]  = utils.fsk2_demod(y_das, bits, Rs, sps, alpha, span, fd);
                        [bits_hat_mvdr, BER_mvdr, pam_rx_mf3, sym_rx3]= utils.fsk2_demod(y_capon, bits, Rs, sps, alpha, span, fd);

                        fprintf('BER: RX SEM INT. %.2f%%  / RX %.2f%% / Capon %.2f%% / DAS %.2f%% \n', ...
                            BER_q*100, BER_in*100, BER_mvdr*100, BER_das*100);

                        % Calcula EVM
                        [EVM_onlynoise,  EVMdB_onlynoise]   = utils.calc_evm_real(sym_rx0,  sym_tx);
                        [EVM_controle,  EVMdB_controle]   = utils.calc_evm_real(sym_rx4,  sym_tx);
                        [EVM_in,  EVMdB_in]                 = utils.calc_evm_real(sym_rx1,  sym_tx);
                        [EVM_das, EVMdB_das]                = utils.calc_evm_real(sym_rx2, sym_tx);
                        [EVM_mvdr,EVMdB_mvdr]               = utils.calc_evm_real(sym_rx3, sym_tx);

                        fprintf('EVM (dB): Rx Sem Int. %.2f | Rx: %.2f | Capon: %.2f | DAS: %.2f\n', ...
                            EVMdB_onlynoise , EVMdB_in, EVMdB_mvdr, EVMdB_das);

                        % Salvando EVM e BER
                        BER(:,iSNR,iISR,iRadius, iCoupling, ii) = [
                            % BER_q;
                            BER_controle;
                            BER_in;
                            BER_mvdr;
                            BER_das
                            ];

                        EVM(:,iSNR,iISR,iRadius, iCoupling, ii) = [
                            % EVMdB_onlynoise;
                            BER_controle;
                            EVMdB_in;
                            EVMdB_mvdr;
                            EVMdB_das
                            ];

                        ii = ii + 1;

                    end
                end
            end
        end
    end
end

BER_mean = mean(BER, 6);   % média ao longo da 6ª dimensão (trials)
EVM_mean = mean(EVM, 6);


%% ======= Gráficos de BER e EVM por SNR, ISR, CF após beamforming =======

methods = ["SINAL ORIGINAL", "SEM BF","MPDR","DAS"];
nMethods = 4;
markers = ["o-","x-","s-","d-","^-","v-","*-","+-"];
colors =  ["red", "green", "blue", "black", "red", "green", "blue", "black"];
line_style =   ["-", "--", ":"]; % linha continua sem acoplamento % linha tracejada com acoplamento

% Criar pasta 'graficos'
outDir = fullfile(pwd, 'graficosBeamforming');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% GRÁFICOS BER

%  GRÁFICOS: Varredura de SNR

for iRadius = 1:nRadius
    for iISR = 1:nISR
        fig = figure('Name', "BER vs SNR | ISR=" + string(range_ISR_dB(iISR)) + " (r = " + string(range_radius(iRadius)) + ")", ...
            'NumberTitle', 'off');
        hold on; grid on;

        legend_strings = strings(0);   % inicializa como array de strings vazio

        for iCoupling = 1:nCoupling 
            for m = 1:nMethods % ignorando o primeiro método
                plot(range_SNR_dB, squeeze(BER_mean(m, :, iISR, iRadius, iCoupling))*100, ...
                    markers(m), 'Color', colors(m),'LineWidth', 2, 'LineStyle', line_style(iCoupling));
            end
        end
        for m = 1:nMethods % ignorando o primeiro método
            legend_strings(end+1) = methods(m);
        end

        xlabel("SNR (dB)",'Interpreter','latex');
        ylabel("BER (%)",'Interpreter','latex');
        ylim([0 60]);
        % title(sprintf('BER vs SNR | ISR=%d dB | Snapshots=%d', ...
        %     range_ISR_dB(iISR), range_snapshots(end)), 'FontSize', 12);

        % legend(legend_strings, 'Location', 'northwest', 'FontSize', 12);
        legend(legend_strings, 'Location', 'best','Interpreter','latex');
        
        fileName = sprintf('BER_vs_SNR_ISR_%ddB_Snap_%d.eps', ...
            range_ISR_dB(iISR), range_snapshots(end));

        exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);
        
        fileName = sprintf('BER_vs_SNR_ISR_%ddB_Snap_%d.tex', ...
            range_ISR_dB(iISR), range_snapshots(end));

        %cleanfigure;
        % Define LaTeX macros for width and height (you will define these in your .tex file)
        matlab2tikz(fullfile(outDir, fileName), 'width', '\figurewidth', 'height', '\figureheight')
    end
end

%  GRÁFICOS: Varredura de ISR

for iRadius = 1:nRadius
    for iSNR = 1:nSNR
        fig = figure('Name', "BER vs ISR | SNR=" + string(range_SNR_dB(iSNR)) + " (r = " + string(range_radius(iRadius)) + ")", ...
            'NumberTitle', 'off');
        hold on; grid on;

        legend_strings = strings(0);   % inicializa como array de strings vazio

        for iCoupling = 1:nCoupling
            for m = 1:nMethods % ignorando o primeiro método
                plot(range_ISR_dB, squeeze(BER_mean(m, iSNR, :, iRadius, iCoupling))*100, ...
                    markers(m), 'Color', colors(m),'LineWidth', 2, 'LineStyle', line_style(iCoupling));
            end
        end
        for m = 1:nMethods % ignorando o primeiro método
            legend_strings(end+1) = methods(m);
        end

        xlabel("ISR (dB)",'Interpreter','latex');
        ylabel("BER (%)",'Interpreter','latex');
        ylim([0 60]);
        % title(sprintf('BER vs SNR | ISR=%d dB | Snapshots=%d', ...
        %     range_ISR_dB(iISR), range_snapshots(end)), 'FontSize', 12);

        % legend(legend_strings, 'Location', 'northwest', 'FontSize', 12);
        legend(legend_strings, 'Location', 'best','Interpreter','latex');
        
        fileName = sprintf('BER_vs_ISR_SNR_%ddB_Snap_%d.eps', ...
            range_SNR_dB(iSNR), range_snapshots(end));

        exportgraphics(fig, fullfile(outDir, fileName), 'Resolution', 300);
        
        fileName = sprintf('BER_vs_ISR_SNR_%ddB_Snap_%d.tex', ...
            range_SNR_dB(iSNR), range_snapshots(end));

        %cleanfigure;
        % Define LaTeX macros for width and height (you will define these in your .tex file)
        matlab2tikz(fullfile(outDir, fileName), 'width', '\figurewidth', 'height', '\figureheight')
    end
end


return;


%% GRÁFICOS EVM

%  GRÁFICOS: Varredura de SNR

for iRadius = 1:nRadius
    for iISR = 1:nISR
        fig = figure('Name', "EVM vs SNR | ISR=" + string(range_ISR_dB(iISR)), ...
            'NumberTitle', 'off');
        hold on; grid on;

        legend_strings = strings(0);   % inicializa como array de strings vazio

        for iRadius = 1:nRadius
            if iRadius == 1
                line_style = '-'; % linha continua sem acoplamento
            else
                line_style = '--'; % linha tracejada com acoplamento
            end
            for m = 2:nMethods
                if m == 2
                    line_style2 = ':';
                else
                    line_style2 = line_style;
                end
                plot(range_SNR_dB, squeeze(EVM_mean(m,:,iISR,iRadius)), markers(m), 'LineWidth', 4, 'LineStyle', line_style(iRadius));
                legend_strings(end+1) = methods(m) + " with CF: " + string(range_radius(iRadius));
            end
        end

        xlabel("SNR (dB)", 'FontSize', 12);
        ylabel("EVM (dB)", 'FontSize', 12);
        title(sprintf('EVM vs SNR | ISR=%d dB | Snapshots=%d', ...
            range_ISR_dB(iISR), range_snapshots(end)), 'FontSize', 12);

        legend(legend_strings, 'Location', 'northwest', 'FontSize', 12);
        ylim([-30 60]);
        fileName = sprintf('EVM_vs_SNR_ISR_%ddB_Snap_%d.eps', ...
            range_ISR_dB(iISR), range_snapshots(end));

        exportgraphics(fig, fullfile(outDir, fileName));
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
