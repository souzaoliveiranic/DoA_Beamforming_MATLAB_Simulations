% =========================================================================
% DoA_KW_UCA_Manifold_Demo
%
% Demonstra a UTILIDADE de usar o MANIFOLD (steering efetivo C a(theta)) em
% vez de ESTIMAR e INVERTER C. O manifold e' ESTIMADO a partir do mesmo
% Chat do self-cal (metodo de compensacao) -> comparacao justa: o mesmo Chat,
% mudando so' INVERTER (D=inv(Chat)) vs USAR COMO MANIFOLD (Chat a(theta)).
%
%   Estrategias (estimador de DoA = Capon/MPDR):
%     - ideal             : sem acoplamento, dicionario a(theta)
%     - acoplado (nao comp): com acoplamento, dicionario IDEAL a(theta)
%     - inversao          : com acoplamento, desacopla D=inv(Chat), dic a(theta)
%     - manifold estimado : com acoplamento, dicionario Chat a(theta)   (sem inverter)
%     - manifold ideal    : referencia com C verdadeiro, dic C a(theta)
%
%   Chat vem do self-cal dirigido por KW (metodo de compensacao robusto).
%
%   Figuras:
%     (1) BEAMPATTERN: 'acoplado' aponta torto; inversao e manifold apontam certo.
%     (2) RMSE de DoA vs SNR: manifold estimado ~ ideal e bate a inversao
%         (que sofre amplificacao/coloracao de ruido), sobretudo em SNR baixa.
% =========================================================================

clear; clc; close all;
warning('off','MATLAB:singularMatrix'); warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:illConditionedMatrix'); warning('off','estimate_C_circulant_uca:smallAlpha');

%% ---- Parametros ----
M = 8; fc = 500e6; c = 3e8; lambda = c/fc; r = 0.1*lambda; theta_sig_deg = 90;
fs = 288000; N = 2100; Rs = 9600; sps = 30; alpha = 0.3; span = 8; fd = 4.8e3;

maxIter = 8; u = 0.5;                 % self-cal (KW) para 'inversao' e 'manifold estimado'
range_SNR_dB = [-9 -6 -3 0 3 6 9 12];
phi_grid_deg = -180:0.5:180;
beta_uca = 2*pi*(0:M-1).'/M;
K = floor(M/2);

% Monte Carlo
n_angles = 24; n_trials = 10; n_real = n_angles*n_trials;
rng(2026,'twister');
gstep = phi_grid_deg(2)-phi_grid_deg(1);
phi_set = round((-180 + 360*rand(1,n_angles))/gstep)*gstep;

outDir = fullfile(pwd,'Manifold Demo Graphs'); if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---- C_true e dicionarios ----
C_true = compute_Ctx_for_R(fc, M, r, 50);
fprintf('C_true (raio=%.2f lambda, cond=%.1f)\n', r/lambda, cond(C_true));
nGrid = numel(phi_grid_deg);
A_id  = zeros(M,nGrid);
for ig = 1:nGrid
    A_id(:,ig) = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_grid_deg(ig));
end
B_cal = C_true*A_id;  B_cal = B_cal ./ vecnorm(B_cal,2,1);   % manifold IDEAL (referencia)

%% =====================================================================
%   PARTE 1: RMSE de DoA vs SNR (Capon, 5 estrategias)
%% =====================================================================
% linhas: 1 ideal | 2 acoplado | 3 inversao | 4 manifold estimado | 5 manifold ideal
nSNR = numel(range_SNR_dB);
sq = zeros(5,nSNR);  cnt = zeros(1,nSNR);
% ganho de array (SNR de saida) na direcao da fonte: 1 ideal | 2 acoplado | 3 inversao | 4 manifold | 5 otimo
sumG = zeros(5,nSNR);  cntG = zeros(5,nSNR);
t0 = tic;
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fprintf('\n== SNR=%+d dB ==\n', SNR_dB);
    for ia = 1:n_angles
        phi_true = phi_set(ia);
        for itr = 1:n_trials
            [~, q, ~, Xsig, ~, Xn] = utils.simulate_fsk_data_uca(M, r, lambda, ...
                phi_true, phi_true+90, theta_sig_deg, theta_sig_deg, ...
                SNR_dB, -100, N, fs, Rs, sps, alpha, span, fd);
            q = q(:);
            X_ideal = Xsig + Xn;  X_coupled = C_true*Xsig + Xn;
            b_hat   = X_coupled*conj(q)/(q'*q);

            R_c = (X_coupled*X_coupled')/N;
            R_i = (X_ideal  *X_ideal'  )/N;

            % --- self-cal (KW) -> Chat (usado por inversao E manifold estimado) ---
            c_fin = selfcal_kw(X_coupled, b_hat, q, M, r, lambda, theta_sig_deg, ...
                               beta_uca, A_id, phi_grid_deg, u, maxIter);
            C_hat = reconstruct_C(c_fin, M);
            D = safe_inv(C_hat);
            Yd = D*X_coupled;  R_d = (Yd*Yd')/N;                 % dado desacoplado
            B_est = C_hat*A_id;  B_est = B_est ./ vecnorm(B_est,2,1);  % manifold ESTIMADO

            phi_ideal = capon_doa(R_i, A_id,  phi_grid_deg);
            phi_acop  = capon_doa(R_c, A_id,  phi_grid_deg);
            phi_inv   = capon_doa(R_d, A_id,  phi_grid_deg);
            phi_mest  = capon_doa(R_c, B_est, phi_grid_deg);
            phi_mid   = capon_doa(R_c, B_cal, phi_grid_deg);

            sq(1,iSNR)=sq(1,iSNR)+wrapTo180(phi_ideal-phi_true)^2;
            sq(2,iSNR)=sq(2,iSNR)+wrapTo180(phi_acop -phi_true)^2;
            sq(3,iSNR)=sq(3,iSNR)+wrapTo180(phi_inv  -phi_true)^2;
            sq(4,iSNR)=sq(4,iSNR)+wrapTo180(phi_mest -phi_true)^2;
            sq(5,iSNR)=sq(5,iSNR)+wrapTo180(phi_mid  -phi_true)^2;

            % --- ganho de array (output SNR) na direcao da fonte ---
            a_n = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_true);  hh = C_true*a_n;
            Gv = [ real(a_n'*a_n); ...                              % ideal a(theta) (ref)
                   abs(a_n'*hh)^2/real(a_n'*a_n); ...               % acoplado  (w=a)
                   abs(a_n'*(D*hh))^2/real(a_n'*(D*D')*a_n); ...    % inversao  (w=D^H a)
                   abs(b_hat'*hh)^2/real(b_hat'*b_hat); ...         % manifold  (w=b_hat)
                   real(hh'*hh) ];                                  % otimo     (w=Ca)
            for kk = 1:5
                if isfinite(Gv(kk)) && Gv(kk) > 0
                    sumG(kk,iSNR)=sumG(kk,iSNR)+Gv(kk); cntG(kk,iSNR)=cntG(kk,iSNR)+1;
                end
            end
            cnt(iSNR)=cnt(iSNR)+1;
        end
    end
    fprintf('  %.0fs\n', toc(t0));
end
rmse = sqrt(sq ./ cnt);
flr = @(x) max(x,1e-3);

labs = {'ideal (sem acopl.)','acoplado (nao comp.)','inversao  D=inv(C_{est})', ...
        'manifold estimado  C_{est}a(\theta)','manifold ideal  C\,a(\theta)'};
fig = figure('Color','w','Position',[80 80 940 580]); hold on; grid on;
% Underlays grossos (referencias) desenhados primeiro; curvas com estilos
% distintos (tracejado/traco-ponto) desenhadas por cima -> mesmo quando
% coincidem, os tracos se intercalam e todas ficam visiveis.
plot(range_SNR_dB, flr(rmse(1,:)), '-',  'Color',[0.55 0.55 0.55],'LineWidth',3.4,'DisplayName',labs{1});                 % ideal (underlay cinza)
plot(range_SNR_dB, flr(rmse(5,:)), '-',  'Color',[0.62 0.85 0.62],'LineWidth',3.4,'DisplayName',labs{5});                 % manifold ideal (underlay verde claro)
plot(range_SNR_dB, flr(rmse(2,:)), 's-', 'Color',[0.75 0.10 0.10],'LineWidth',1.8,'MarkerFaceColor',[0.75 0.10 0.10],'MarkerSize',7,'DisplayName',labs{2});  % acoplado
plot(range_SNR_dB, flr(rmse(3,:)), 'd--','Color',[0.10 0.25 0.85],'LineWidth',1.8,'MarkerFaceColor',[0.10 0.25 0.85],'MarkerSize',7,'DisplayName',labs{3});  % inversao (tracejada)
plot(range_SNR_dB, flr(rmse(4,:)), 'o-.','Color',[0.00 0.45 0.10],'LineWidth',1.8,'MarkerFaceColor',[0.00 0.45 0.10],'MarkerSize',7,'DisplayName',labs{4});  % manifold estimado (traco-ponto)
set(gca,'YScale','log'); xlabel('SNR (dB)'); ylabel('RMSE de DoA (graus)');
xticks(range_SNR_dB); title(sprintf(['DoA (Capon): inversao \\equiv manifold (equivalentes p/ DoA); ambos corrigem o vies\n' ...
    'raio=%.2f\\lambda, %d azimutes x %d realizacoes,  C_{est}: self-cal KW'], r/lambda, n_angles, n_trials));
legend('Location','northeast');
exportgraphics(fig, fullfile(outDir,'manifold_doa_rmse_vs_snr.png'),'Resolution',180);

%% =====================================================================
%   PARTE 2: GANHO DE ARRAY (SNR de saida) vs SNR de entrada
%% =====================================================================
% A diferenca manifold vs inversao NAO aparece na DoA (invariante a
% transformacao invertivel). Aparece aqui: inverter (D=inv(Chat)) aplicado
% AOS DADOS amplifica o ruido; o manifold (filtro casado w=b_hat) nao.
% Ganho de array = |w^H h|^2 / ||w||^2 (h = C a, resposta real da fonte).
GdB = 10*log10( sumG ./ max(cntG,1) );   % 5 x nSNR
labG = {'ideal  a(\theta)  [ref 0 dB]','acoplado  (w=a)','inversao  (w=D^Ha)', ...
        'manifold  (w=b_{hat})','otimo  (w=Ca)'};
fig = figure('Color','w','Position',[80 80 940 580]); hold on; grid on;
plot(range_SNR_dB, GdB(5,:), '-',  'Color',[0.2 0.2 0.2],'LineWidth',3.0,'DisplayName',labG{5});                       % otimo (teto)
plot(range_SNR_dB, GdB(1,:), '--', 'Color',[0.55 0.55 0.55],'LineWidth',1.6,'DisplayName',labG{1});                    % ideal (ref)
plot(range_SNR_dB, GdB(2,:), 's-', 'Color',[0.75 0.10 0.10],'LineWidth',1.8,'MarkerFaceColor',[0.75 0.10 0.10],'MarkerSize',7,'DisplayName',labG{2});  % acoplado
plot(range_SNR_dB, GdB(3,:), 'd--','Color',[0.10 0.25 0.85],'LineWidth',2.0,'MarkerFaceColor',[0.10 0.25 0.85],'MarkerSize',7,'DisplayName',labG{3});  % inversao
plot(range_SNR_dB, GdB(4,:), 'o-', 'Color',[0.00 0.45 0.10],'LineWidth',2.4,'MarkerFaceColor',[0.00 0.45 0.10],'MarkerSize',8,'DisplayName',labG{4});  % manifold
xlabel('SNR de entrada (dB)'); ylabel('ganho de array  10log_{10}(|w^Hh|^2/||w||^2)  (dB)');
xticks(range_SNR_dB); title(sprintf(['Beamforming: ganho de array vs SNR  (raio=%.2f\\lambda)\n' ...
    'manifold (casado, sem inverter)  \\geq  inversao (amplifica ruido)'], r/lambda));
legend('Location','best');
exportgraphics(fig, fullfile(outDir,'manifold_ganho_array_vs_snr.png'),'Resolution',180);

%% ---- Resumo ----
fprintf('\n===== RMSE de DoA (graus) =====\n');
fprintf('%-26s', 'SNR(dB)'); fprintf('%8d', range_SNR_dB); fprintf('\n');
for k=1:5, fprintf('%-26s', labs{k}); fprintf('%8.3f', rmse(k,:)); fprintf('\n'); end
fprintf('\nFiguras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function phi = capon_doa(R, dict, phi_grid)
    M = size(R,1);
    R = R + 1e-6*trace(R)/M*eye(M);
    Rinv = R\eye(M);
    den = real(sum(conj(dict).*(Rinv*dict), 1));   % P = 1/den ; pico = min(den)
    [~,ip] = min(max(den,eps));
    phi = phi_grid(ip);
end

function c_fin = selfcal_kw(X_coupled, b_hat_orig, q, M, r, lambda, theta_deg, ...
                            beta_uca, A_dict, phi_grid_deg, u, maxIter)
    phi0 = doa_estimate(X_coupled, 'KW', q, r, lambda, beta_uca, A_dict, phi_grid_deg);
    C_hat = estimate_C_circulant_uca(b_hat_orig, utils.steering_vec_uca(M,r,lambda,theta_deg,phi0), M);
    for it = 1:maxIter
        D = safe_inv(C_hat);  Y = D*X_coupled;
        phi = doa_estimate(Y, 'KW', q, r, lambda, beta_uca, A_dict, phi_grid_deg);
        a_hat = utils.steering_vec_uca(M, r, lambda, theta_deg, phi);
        C_ls = estimate_C_circulant_uca(b_hat_orig, a_hat, M);
        C_hat = (1-u)*C_hat + u*C_ls;
    end
    c_fin = C_hat(1, 2:floor(M/2)+1).';
end

function C = reconstruct_C(c, M)
    K = floor(M/2); is_even = (mod(M,2)==0); c = c(:);
    if is_even, first_row = [1, c(1:K-1).', c(K), flip(c(1:K-1).')];
    else,       first_row = [1, c(1:K).', flip(c(1:K).')]; end
    C = zeros(M);
    for i = 1:M, C(i,:) = circshift(first_row, [0, i-1]); end
end

function D = safe_inv(C)
    M = size(C,1);
    ws = warning('off','MATLAB:singularMatrix'); wn = warning('off','MATLAB:nearlySingularMatrix');
    cleanupObj = onCleanup(@() warning([ws wn]));
    rc = rcond(C);
    if isfinite(rc) && rc > 1e-12, D = inv(C); if all(isfinite(D(:))), return; end, end
    mu = 1e-6*(norm(C,'fro')^2/M + eps);  D = (C'*C + mu*eye(M))\C';
    if ~all(isfinite(D(:))), D = eye(M); end
end

function phi_hat = doa_estimate(X, method, q, r, lambda, beta_uca, A_dict, phi_grid_deg)
    M = size(X,1);
    switch upper(method)
        case 'KW'
            [~, phi] = doa_kw_uca(X, q(:).', r, lambda, beta_uca);  phi_hat = phi(1);
        case 'CAPON'
            R = (X*X')/size(X,2); R = R + 1e-6*trace(R)/M*eye(M); Rinv = R\eye(M);
            den = real(sum(conj(A_dict).*(Rinv*A_dict),1)); [~,ip]=max(1./max(den,eps)); phi_hat=phi_grid_deg(ip);
        otherwise
            error('doa_estimate: metodo %s', method);
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
