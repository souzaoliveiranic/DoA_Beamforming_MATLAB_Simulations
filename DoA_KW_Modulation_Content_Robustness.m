% =========================================================================
% DoA_KW_Modulation_Content_Robustness.m
%
% Estuda a robustez do estimador KW (known-waveform) de DoA em UCA compacta
% com relacao a DOIS fatores ligados ao sinal de referencia (SOI):
%
%   (A) MODULACAO do SOI:
%         FSK2, FSK4, QPSK, QAM16, QAM64, e um sinal ALEATORIO (gaussiano,
%         sem estrutura) usado como referencia teorica (upper bound).
%
%   (B) ESTRUTURA DE CONTEUDO do preambulo:
%         - 'random'   : sequencia de simbolos totalmente aleatoria (baseline)
%         - 'repeated' : uma sequencia curta de comprimento L_seq repetida
%                        ate preencher os K snapshots (preambulos reais).
%         Varia-se L_seq para mapear o trade-off entre comprimento total de
%         observacao (K) e diversidade temporal efetiva (K/L_seq).
%
% Metrica: RMSE do erro angular ESFERICO (geodesico) 2D, como no artigo
% principal. O interferente e' tratado como perturbacao (FM narrowband por
% padrao). O foco esta no estimador KW; os classicos (DAS/Capon/MUSIC) sao
% calculados como referencia.
%
% Este script NAO altera utils.m. Reaproveita:
%   utils.gen_linmod_interferer, utils.steering_vec_uca,
%   utils.element_delays_uca, utils.gen_ce_nb_noise, utils.compute_Ctx_for_R
%   doa_kw_uca (funcao externa), spherical_angular_distance (local, no fim)
%
% =========================================================================
clear; clc; close all;

%% ----------------------------- Array -----------------------------------
M      = 8;            % nº de elementos do UCA
fc     = 500e6;        % Hz
c      = 3e8;
lambda = c/fc;

Z0     = 50;           % impedancia de referencia
beta   = 2*pi*(0:M-1)'/M;

%% --------------------------- Sinal / banda -----------------------------
fs     = 288000;       % taxa de amostragem (Hz)
N      = 21000;        % nº de amostras geradas
N_DOA  = 2000;         % snapshots usados na estimacao (sample support K)

Rs     = 9600;         % taxa de simbolos [sym/s]
sps    = 30;           % amostras por simbolo
alpha  = 0.3;          % roll-off do RRC
span   = 8;            % span do RRC (simbolos)
fd     = 4.8e3;        % desvio de frequencia (FSK) [Hz]

%% ------------------------- Eixos do experimento ------------------------
% (A) Modulacoes do SOI a avaliar
mod_list = ["FSK2","FSK4","QPSK","QAM16","QAM64","RANDOM"];
nMod = numel(mod_list);

% (B) Estrutura de conteudo do preambulo
%   content_mode = 'random'   -> sequencia aleatoria (L_seq ignorado)
%   content_mode = 'repeated' -> sequencia de L_seq simbolos repetida
content_mode = 'repeated';     % 'random' | 'repeated'
% Periodos de repeticao a varrer (em SIMBOLOS) no modo 'repeated'.
% L_seq grande -> proximo de aleatorio; L_seq pequeno -> muita repeticao.
range_Lseq = [1 2 4 8 16 32 64 128];
if strcmpi(content_mode,'random')
    range_Lseq = NaN;   % um unico ponto, sem repeticao
end
nLseq = numel(range_Lseq);

% Robustez: varre SNR (interferente e coupling fixos para isolar o efeito
% da modulacao / conteudo). Ajuste conforme a figura desejada.
range_SNR_dB   = -6:3:6;
ISR_dB_fixed   = -40;            % interferente com mesma potencia do SOI
range_coupling = [0 1];        % 0 = sem coupling, 1 = com coupling
radius_factor  = 0.22;         % raio do UCA em unidades de lambda
radius         = radius_factor * lambda;

interf_type    = 'fm';         % interferente: 'fm' (narrowband CE)

% Monte Carlo
n_rand_angles  = 200;          % pares angulares (SOI, interferente) por ponto
min_sep_deg    = 10;           % separacao esferica minima (graus)
theta_min_rand = 10;
theta_max_rand = 90;

nSNR      = numel(range_SNR_dB);
nCoupling = numel(range_coupling);
nMethods  = 4;
methods   = ["KW","DAS","Capon","MUSIC"];

%% --------------------- Pares angulares (fixos) -------------------------
rng(12345,'twister');
sig_dirs = zeros(n_rand_angles, 2);   % [phi, theta]
int_dirs = zeros(n_rand_angles, 2);
for kk = 1:n_rand_angles
    ps = -180 + 360*rand;
    ts = theta_min_rand + (theta_max_rand - theta_min_rand)*rand;
    while true
        pi_ = -180 + 360*rand;
        ti_ = theta_min_rand + (theta_max_rand - theta_min_rand)*rand;
        if spherical_angular_distance(ts, ps, ti_, pi_) >= min_sep_deg
            break;
        end
    end
    sig_dirs(kk,:) = [ps, ts];
    int_dirs(kk,:) = [pi_, ti_];
end

%% --------------------- Grade de busca 2D (classicos) -------------------
theta_scan_grid = 1:1:90;        % evita zenith
phi_scan_grid   = -180:0.5:180;
[PhiScanG, ThetaScanG] = meshgrid(phi_scan_grid, theta_scan_grid);
theta_scan_flat = ThetaScanG(:);
phi_scan_flat   = PhiScanG(:);
nGridPts        = numel(theta_scan_flat);

fprintf('Pre-computando steering vectors da grade 2D (%d pontos)...\n', nGridPts);
tic_pre = tic;
A_scan = zeros(M, nGridPts);
for kk = 1:nGridPts
    A_scan(:,kk) = utils.steering_vec_uca(M, radius, lambda, ...
                                          theta_scan_flat(kk), phi_scan_flat(kk));
end
fprintf('  ...concluido em %.2f s\n', toc(tic_pre));

%% --------------------- Matriz de acoplamento ---------------------------
Ctx = compute_Ctx_for_R(fc, M, radius, Z0);
Coupling_matrix_on = Ctx;     % aplicada quando coupling==1

%% --------------------- Armazenamento -----------------------------------
% RMSE(metodo, SNR, coupling, modulacao, Lseq, par_angular)
RMSE = zeros(nMethods, nSNR, nCoupling, nMod, nLseq, n_rand_angles);

TotalSim = nSNR * nCoupling * nMod * nLseq * n_rand_angles;
iTotal = 0;
sweep_t0 = tic;
fprintf('\n=== INICIO DO LOOP PRINCIPAL ===\n');
fprintf('Total de simulacoes: %d (modo de conteudo: %s)\n', TotalSim, content_mode);
fprintf('================================\n\n');

%% --------------------- Loop principal ----------------------------------
for iMod = 1:nMod
    modtype = mod_list(iMod);

    for iLseq = 1:nLseq
        Lseq = range_Lseq(iLseq);

        for iSNR = 1:nSNR
            SNR_dB = range_SNR_dB(iSNR);

            for iCoupling = 1:nCoupling
                Coupling = range_coupling(iCoupling);
                if Coupling == 0
                    Cmat = eye(M);
                else
                    Cmat = Coupling_matrix_on;
                end

                for iPair = 1:n_rand_angles
                    iTotal = iTotal + 1;

                    phi_sig = sig_dirs(iPair,1);  theta_sig = sig_dirs(iPair,2);
                    phi_int = int_dirs(iPair,1);  theta_int = int_dirs(iPair,2);

                    % --- progresso ---
                    if mod(iTotal, 200) == 0 || iTotal == 1
                        elapsed = toc(sweep_t0);
                        if iTotal > 1
                            eta_sec = elapsed*(TotalSim-iTotal)/(iTotal-1);
                        else
                            eta_sec = 0;
                        end
                        fprintf(['Sim %d/%d (%.1f%%) | mod=%s Lseq=%g SNR=%+d Coup=%d | ' ...
                                 'elapsed %s ETA %s\n'], ...
                                iTotal, TotalSim, 100*iTotal/TotalSim, ...
                                modtype, Lseq, SNR_dB, Coupling, ...
                                fmt_hms(elapsed), fmt_hms(eta_sec));
                    end

                    % ============================================================
                    % 1) Gera SOI baseband com modulacao + estrutura de conteudo
                    % ============================================================
                    q = gen_soi_waveform(modtype, content_mode, Lseq, N, Rs, sps, alpha, span, fd, fs);

                    % 2) Interferente (FM narrowband CE)
                    v = utils.gen_ce_nb_noise(N, fs, 2*fd, 0);

                    % 3) Normaliza potencias conforme SNR/ISR
                    sigma_s2 = 1;
                    sigma_i2 = sigma_s2 * 10^(ISR_dB_fixed/10);
                    sigma_n2 = sigma_s2 / 10^(SNR_dB/10);

                    q = q ./ sqrt(mean(abs(q).^2) + eps);  q = sqrt(sigma_s2)*q;
                    v = v ./ sqrt(mean(abs(v).^2) + eps);  v = sqrt(sigma_i2)*v;

                    % 4) Monta sinais no array (atraso + fase espacial)
                    taus_sig = utils.element_delays_uca(M, radius, theta_sig, phi_sig, c);
                    taus_int = utils.element_delays_uca(M, radius, theta_int, phi_int, c);

                    Xsig = zeros(M,N); Xint = zeros(M,N);
                    for m = 1:M
                        q_del = delayseq(q, taus_sig(m), fs);
                        v_del = delayseq(v, taus_int(m), fs);
                        Xsig(m,:) = q_del(:).' * exp(-1j*2*pi*fc*taus_sig(m));
                        Xint(m,:) = v_del(:).' * exp(-1j*2*pi*fc*taus_int(m));
                    end
                    Xn = sqrt(sigma_n2/2)*(randn(M,N)+1j*randn(M,N));

                    % 5) Aplica coupling e soma
                    X = Cmat*Xsig + Cmat*Xint + Xn;

                    K = N_DOA;

                    % ============================================================
                    % Estimadores
                    % ============================================================
                    % --- KW (referencia = SOI conhecido q) ---
                    [theta_KW, phi_KW] = doa_kw_uca(X(:,1:K), q(1:K).', radius, lambda, beta);

                    % --- Classicos (busca 2D vetorizada) ---
                    Rxx = (X(:,1:K)*X(:,1:K)')/K;
                    delta = 1e-3*trace(Rxx)/M;
                    Rxx_dl = Rxx + delta*eye(M);
                    Rinv = inv(Rxx_dl);
                    [eigvec, eigval] = eig(Rxx_dl);
                    [~, idx] = sort(diag(eigval),'descend');
                    En = eigvec(:, idx(3:end));   % Ksrc=2
                    EnEnH = En*En';

                    P_DAS   = abs(  sum(conj(A_scan).*(Rxx   *A_scan),1) );
                    P_MVDR  = 1./max(real(sum(conj(A_scan).*(Rinv *A_scan),1)),eps);
                    P_MUSIC = 1./max(real(sum(conj(A_scan).*(EnEnH*A_scan),1)),eps);

                    [~,i_das  ] = max(P_DAS);
                    [~,i_mvdr ] = max(P_MVDR);
                    [~,i_music] = max(P_MUSIC);

                    % ============================================================
                    % Erro esferico 2D
                    % ============================================================
                    e_KW    = spherical_angular_distance(theta_sig,phi_sig,theta_KW,            phi_KW);
                    e_DAS   = spherical_angular_distance(theta_sig,phi_sig,theta_scan_flat(i_das),  phi_scan_flat(i_das));
                    e_MVDR  = spherical_angular_distance(theta_sig,phi_sig,theta_scan_flat(i_mvdr), phi_scan_flat(i_mvdr));
                    e_MUSIC = spherical_angular_distance(theta_sig,phi_sig,theta_scan_flat(i_music),phi_scan_flat(i_music));

                    RMSE(:,iSNR,iCoupling,iMod,iLseq,iPair) = [e_KW; e_DAS; e_MVDR; e_MUSIC];
                end % iPair
            end % iCoupling
        end % iSNR
    end % iLseq
end % iMod

fprintf('\nLoop concluido em %s\n', fmt_hms(toc(sweep_t0)));

% RMSE quadratico medio sobre os pares angulares (ultima dim)
RMSE_mean = sqrt(mean(RMSE.^2, 6));   % (metodo, SNR, coupling, mod, Lseq)

%% =======================================================================
%  GRAFICOS
%  Dois cortes principais:
%   (1) RMSE do KW vs SNR, uma curva por MODULACAO (Lseq fixo / random)
%   (2) RMSE do KW vs Lseq (diversidade temporal), uma curva por MODULACAO
% =======================================================================
colors = lines(nMod);

% ---- Figura 1: KW vs SNR, por modulacao (sem coupling, Lseq de referencia) ----
iLseq_ref = nLseq;   % maior Lseq (mais proximo de aleatorio) como referencia
figure('Name','KW vs SNR por modulacao'); hold on; grid on;
for iMod = 1:nMod
    y = squeeze(RMSE_mean(1, :, 1, iMod, iLseq_ref));  % metodo=1 (KW), coupling=0
    plot(range_SNR_dB, y, 'o-', 'Color', colors(iMod,:), 'LineWidth', 1.4, ...
        'DisplayName', mod_list(iMod));
end
set(gca,'YScale','log'); ylim([0.5 100]);
xlabel('SNR (dB)'); ylabel('RMSE (degrees)');
title(sprintf('KW DoA — RMSE vs SNR by modulation (no coupling, L_{seq}=%g)', range_Lseq(iLseq_ref)));
legend('Location','best');

% ---- Figura 2: KW vs Lseq, por modulacao (so no modo 'repeated') ----
if strcmpi(content_mode,'repeated')
    iSNR_ref = find(range_SNR_dB==0, 1); if isempty(iSNR_ref), iSNR_ref = ceil(nSNR/2); end
    figure('Name','KW vs Lseq por modulacao'); hold on; grid on;
    for iMod = 1:nMod
        y = squeeze(RMSE_mean(1, iSNR_ref, 1, iMod, :));  % KW, coupling=0, SNR ref
        plot(range_Lseq, y, 'o-', 'Color', colors(iMod,:), 'LineWidth', 1.4, ...
            'DisplayName', mod_list(iMod));
    end
    set(gca,'XScale','log','YScale','log'); ylim([0.5 100]);
    xlabel('Repetition period L_{seq} (symbols)'); ylabel('RMSE (degrees)');
    title(sprintf('KW DoA — RMSE vs preamble repetition period (SNR=%d dB, no coupling)', range_SNR_dB(iSNR_ref)));
    legend('Location','best');
end

% ---- Figura 3: comparacao KW vs classicos para uma modulacao (sanity) ----
iMod_ref = find(mod_list=="FSK2",1); if isempty(iMod_ref), iMod_ref=1; end
figure('Name','KW vs classicos (FSK2)'); hold on; grid on;
ls = {'-','--'};
for iCoupling = 1:nCoupling
    for mth = 1:nMethods
        y = squeeze(RMSE_mean(mth, :, iCoupling, iMod_ref, iLseq_ref));
        plot(range_SNR_dB, y, ['o' ls{iCoupling}], 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s (coup=%d)', methods(mth), range_coupling(iCoupling)));
    end
end
set(gca,'YScale','log'); ylim([0.5 100]);
xlabel('SNR (dB)'); ylabel('RMSE (degrees)');
title(sprintf('All methods — %s', mod_list(iMod_ref)));
legend('Location','best');

%% ---- Tabela resumo no console ----
fprintf('\n===== RESUMO: RMSE do KW (graus), sem coupling, Lseq=%g =====\n', range_Lseq(iLseq_ref));
fprintf('%-8s', 'Mod\\SNR');
for iSNR=1:nSNR, fprintf('%8d', range_SNR_dB(iSNR)); end
fprintf('\n');
for iMod=1:nMod
    fprintf('%-8s', mod_list(iMod));
    for iSNR=1:nSNR
        fprintf('%8.2f', RMSE_mean(1,iSNR,1,iMod,iLseq_ref));
    end
    fprintf('\n');
end

%% =======================================================================
%  FUNCOES LOCAIS
% =======================================================================
function q = gen_soi_waveform(modtype, content_mode, Lseq, N, Rs, sps, alpha, span, fd, fs)
    % GEN_SOI_WAVEFORM  Gera o SOI baseband (N x 1) com modulacao e estrutura
    % de conteudo parametrizaveis.
    %
    %   modtype: "FSK2"|"FSK4"|"QPSK"|"QAM16"|"QAM64"|"RANDOM"
    %   content_mode: 'random' | 'repeated'
    %   Lseq: periodo de repeticao em simbolos (usado se 'repeated')
    %
    % Estrategia de conteudo:
    %   - Gera a sequencia de SIMBOLOS conforme a modulacao.
    %   - Se 'repeated', constroi um bloco-base de Lseq simbolos e o repete
    %     ate cobrir o numero de simbolos necessario, ANTES do pulse-shaping.

    Nsym = ceil(N/sps) + 2*span + 4;

    modtype = upper(string(modtype));

    % ---- 1) Gera a sequencia de simbolos base (complexa) ----
    switch modtype
        case "FSK2"
            % FSK e' modulacao de frequencia: tratamos no dominio de simbolos
            % como bits +-1 e modulamos FM apos repetir/expandir.
            sym = gen_symbols_pam(2, Nsym);
        case "FSK4"
            sym = gen_symbols_pam(4, Nsym);   % 4 niveis: -3,-1,+1,+3 (normalizado)
        case "QPSK"
            sym = gen_symbols_psk(4, Nsym);
        case "QAM16"
            sym = gen_symbols_qam(16, Nsym);
        case "QAM64"
            sym = gen_symbols_qam(64, Nsym);
        case "RANDOM"
            % Sinal gaussiano complexo IID por amostra (sem estrutura de simbolo).
            % Caso ideal: referencia teorica.
            q = (randn(N,1)+1j*randn(N,1))/sqrt(2);
            % conteudo 'repeated' tambem se aplica: repete um bloco de Lseq*sps amostras
            if strcmpi(content_mode,'repeated') && ~isnan(Lseq)
                block = Lseq*sps;
                base  = (randn(block,1)+1j*randn(block,1))/sqrt(2);
                nrep  = ceil(N/block);
                q = repmat(base, nrep, 1);
                q = q(1:N);
            end
            return;
        otherwise
            error('Modulacao nao suportada: %s', modtype);
    end

    % ---- 2) Aplica estrutura de conteudo (repeticao) no dominio de simbolos ----
    if strcmpi(content_mode,'repeated') && ~isnan(Lseq)
        base = sym(1:min(Lseq,numel(sym)));
        nrep = ceil(Nsym/numel(base));
        sym  = repmat(base, nrep, 1);
        sym  = sym(1:Nsym);
    end

    % ---- 3) Pulse-shaping / modulacao para banda base ----
    switch modtype
        case {"FSK2","FSK4"}
            % FM complexa: integra a sequencia PAM (mesma logica de fsk2_mod)
            imp = upsample(real(sym), sps);
            rrc = rcosdesign(alpha, span, sps, 'sqrt');
            m   = filter(rrc, 1, imp);
            m   = m ./ (max(abs(m)) + eps);
            phi = 2*pi*cumsum(fd*m)/fs;
            s   = exp(1j*phi);
        otherwise
            % PSK / QAM: pulse-shaping RRC linear
            imp = upsample(sym, sps);
            rrc = rcosdesign(alpha, span, sps, 'sqrt');
            s   = filter(rrc, 1, imp);
    end

    % ---- 4) Ajusta comprimento para N ----
    gd = span*sps/2;
    if numel(s) > gd, s = s(gd+1:end); end
    if numel(s) < N, s(end+1:N) = 0; else, s = s(1:N); end
    q = s(:);

    % ---- 5) Normaliza potencia para ~1 ----
    q = q ./ sqrt(mean(abs(q).^2) + eps);
end

function sym = gen_symbols_pam(Mord, Nsym)
    % PAM real centrado e normalizado (para FSK)
    levels = -(Mord-1):2:(Mord-1);          % e.g. M=4 -> [-3 -1 1 3]
    idx = randi(Mord, Nsym, 1);
    sym = levels(idx).';
    sym = sym ./ sqrt(mean(abs(levels).^2)); % potencia unitaria
    sym = sym + 0j;
end

function sym = gen_symbols_psk(Mord, Nsym)
    idx = randi(Mord, Nsym, 1) - 1;
    sym = exp(1j*(2*pi*idx/Mord + pi/Mord));  % PSK com offset
end

function sym = gen_symbols_qam(Mord, Nsym)
    data = randi(Mord, Nsym, 1) - 1;
    sym = qammod(data, Mord, 'gray', 'UnitAveragePower', true);
    sym = sym(:);
end

function s = fmt_hms(sec)
    if ~isfinite(sec) || sec < 0, s = '--:--:--'; return; end
    h = floor(sec/3600); m = floor((sec-3600*h)/60); s2 = floor(sec-3600*h-60*m);
    s = sprintf('%02d:%02d:%02d', h, m, s2);
end

function d_deg = spherical_angular_distance(theta1_deg, phi1_deg, theta2_deg, phi2_deg)
    % Distancia angular esferica (geodesica) entre (theta,phi) em graus.
    t1=deg2rad(theta1_deg); p1=deg2rad(phi1_deg);
    t2=deg2rad(theta2_deg); p2=deg2rad(phi2_deg);
    u1=[sin(t1)*cos(p1); sin(t1)*sin(p1); cos(t1)];
    u2=[sin(t2)*cos(p2); sin(t2)*sin(p2); cos(t2)];
    cval = max(min(u1.'*u2, 1), -1);
    d_deg = rad2deg(acos(cval));
end

function Ctx = compute_Ctx_for_R(fc, M, R, Z0)
    % Monta UCA e calcula Z via S-parameters, depois Ctx.
    c = 3e8; lambda = c/fc;
    mp = dipole; mp.Length = 0.5*lambda; mp.Width = 0.01*lambda;
    uca = circularArray; uca.Element = mp; uca.NumElements = M; uca.Radius = R;
    Sobj = sparameters(uca, fc);
    S_matrix = Sobj.Parameters(:,:,1);
    Z_matrix = s2z(S_matrix, Z0);
    Zg = Z0; Zself = diag(Z_matrix); denom = Zself + Zg;
    Ctx = eye(M);
    for j = 1:M
        for i = 1:M
            if i ~= j, Ctx(i,j) = Z_matrix(i,j)/denom(j); end
        end
    end
end
