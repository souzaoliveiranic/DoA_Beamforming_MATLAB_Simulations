% =========================================================================
% DoA_KW_UCA_Manifold_Demo_PorMetodo
%
% Versao POR METODO da demonstracao manifold vs inverter C. Cada metodo de
% DoA (KW/DAS/Capon/MUSIC) dirige seu PROPRIO self-cal (init KW, passo u) ->
% Chat_m. Um painel por metodo (2x2).
%
%   Figura 1 (DoA RMSE vs SNR), por metodo:
%     - ideal        : metodo no dado SEM acoplamento
%     - nao compensado: metodo no dado ACOPLADO
%     - compensado   : metodo no dado DESACOPLADO (D_m = inv(Chat_m))
%
%   Figura 2 (ganho de array vs SNR), por metodo:
%     - otimo (w=Ca) e manifold (w=b_hat): referencias (independem do metodo)
%     - acoplado (w=a): referencia
%     - inversao (w=D_m^H a): usa o Chat do self-cal DESTE metodo
%   Mostra, por metodo, que inverter Chat_m fica ABAIXO do manifold (ruido).
% =========================================================================

clear; clc; close all;
warning('off','MATLAB:singularMatrix'); warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:illConditionedMatrix'); warning('off','estimate_C_circulant_uca:smallAlpha');

%% ---- Parametros ----
M = 8; fc = 500e6; c = 3e8; lambda = c/fc; r = 0.1*lambda; theta_sig_deg = 90;
fs = 288000; N = 2100; Rs = 9600; sps = 30; alpha = 0.3; span = 8; fd = 4.8e3;

maxIter = 8; u = 0.5;
range_SNR_dB = [-9 -6 -3 0 3 6 9 12];
methods = {'KW','DAS','CAPON','MUSIC'};  nM = numel(methods);
phi_grid_deg = -180:0.5:180;
beta_uca = 2*pi*(0:M-1).'/M;
K = floor(M/2);

% Monte Carlo
n_angles = 20; n_trials = 6; n_real = n_angles*n_trials;
rng(2026,'twister');
gstep = phi_grid_deg(2)-phi_grid_deg(1);
phi_set = round((-180 + 360*rand(1,n_angles))/gstep)*gstep;

outDir = fullfile(pwd,'Manifold Demo PorMetodo Graphs'); if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---- C_true e dicionario ----
C_true = compute_Ctx_for_R(fc, M, r, 50);
fprintf('C_true (raio=%.2f lambda, cond=%.1f)\n', r/lambda, cond(C_true));
nGrid = numel(phi_grid_deg);
A_id = zeros(M,nGrid);
for ig = 1:nGrid, A_id(:,ig) = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_grid_deg(ig)); end

%% ---- Acumuladores ----
nSNR = numel(range_SNR_dB);
% DoA: sq(metodo, estrategia{1 ideal,2 naocomp,3 comp}, SNR)
sq = zeros(nM,3,nSNR);
% Beamforming: inversao por metodo + referencias (independem do metodo)
sumG_inv = zeros(nM,nSNR);  cntG_inv = zeros(nM,nSNR);
sumG_ideal=zeros(1,nSNR); sumG_acop=zeros(1,nSNR); sumG_man=zeros(1,nSNR); sumG_opt=zeros(1,nSNR);
cnt = zeros(1,nSNR);

t0 = tic;
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fprintf('\n== SNR=%+d dB ==\n', SNR_dB);
    for ia = 1:n_angles
        phi_true = phi_set(ia);
        a = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_true);  h = C_true*a;
        for itr = 1:n_trials
            [~, q, ~, Xsig, ~, Xn] = utils.simulate_fsk_data_uca(M, r, lambda, ...
                phi_true, phi_true+90, theta_sig_deg, theta_sig_deg, ...
                SNR_dB, -100, N, fs, Rs, sps, alpha, span, fd);
            q = q(:);
            X_ideal = Xsig + Xn;  X_coupled = C_true*Xsig + Xn;
            b_hat = X_coupled*conj(q)/(q'*q);

            % --- referencias de beamforming (independem do metodo) ---
            sumG_ideal(iSNR)=sumG_ideal(iSNR)+real(a'*a);
            sumG_acop(iSNR) =sumG_acop(iSNR) +abs(a'*h)^2/real(a'*a);
            sumG_man(iSNR)  =sumG_man(iSNR)  +abs(b_hat'*h)^2/max(real(b_hat'*b_hat),eps);
            sumG_opt(iSNR)  =sumG_opt(iSNR)  +real(h'*h);

            % --- por metodo: self-cal proprio -> Chat_m ---
            for m = 1:nM
                c_m = selfcal_drive(methods{m}, X_coupled, b_hat, q, M, r, lambda, ...
                                    theta_sig_deg, beta_uca, A_id, phi_grid_deg, u, maxIter);
                D_m = safe_inv(reconstruct_C(c_m, M));

                % DoA: ideal / nao-comp / compensado (dado desacoplado)
                phi_id = doa_estimate(X_ideal,       methods{m}, q, r, lambda, beta_uca, A_id, phi_grid_deg);
                phi_ac = doa_estimate(X_coupled,     methods{m}, q, r, lambda, beta_uca, A_id, phi_grid_deg);
                phi_cp = doa_estimate(D_m*X_coupled, methods{m}, q, r, lambda, beta_uca, A_id, phi_grid_deg);
                sq(m,1,iSNR)=sq(m,1,iSNR)+wrapTo180(phi_id-phi_true)^2;
                sq(m,2,iSNR)=sq(m,2,iSNR)+wrapTo180(phi_ac-phi_true)^2;
                sq(m,3,iSNR)=sq(m,3,iSNR)+wrapTo180(phi_cp-phi_true)^2;

                % Beamforming: inversao com o Chat deste metodo
                Ginv = abs(a'*(D_m*h))^2/real(a'*(D_m*D_m')*a);
                if isfinite(Ginv) && Ginv>0
                    sumG_inv(m,iSNR)=sumG_inv(m,iSNR)+Ginv; cntG_inv(m,iSNR)=cntG_inv(m,iSNR)+1;
                end
            end
            cnt(iSNR)=cnt(iSNR)+1;
        end
    end
    fprintf('  %.0fs\n', toc(t0));
end

%% ---- Consolidacao ----
rmse = sqrt(sq ./ reshape(cnt,1,1,nSNR));      % nM x 3 x nSNR
flr = @(x) max(x,1e-3);
G_inv  = 10*log10(sumG_inv ./ max(cntG_inv,1));     % nM x nSNR
G_ideal= 10*log10(sumG_ideal./cnt);
G_acop = 10*log10(sumG_acop ./cnt);
G_man  = 10*log10(sumG_man  ./cnt);
G_opt  = 10*log10(sumG_opt  ./cnt);

%% =================== FIGURA 1: DoA RMSE por metodo ===================
% ylim comum aos 4 paineis
vv = flr(rmse(:)); ylD = [min(vv)*0.8, max(vv)*1.3];
fig = figure('Color','w','Position',[60 60 1200 820]);
for m = 1:nM
    subplot(2,2,m); hold on; grid on;
    plot(range_SNR_dB, flr(squeeze(rmse(m,1,:))), '--','Color',[0.35 0.35 0.35],'LineWidth',1.8,'DisplayName','ideal (sem acopl.)');
    plot(range_SNR_dB, flr(squeeze(rmse(m,2,:))), 's-','Color',[0.75 0.10 0.10],'LineWidth',1.8,'MarkerFaceColor',[0.75 0.10 0.10],'MarkerSize',6,'DisplayName','nao compensado');
    plot(range_SNR_dB, flr(squeeze(rmse(m,3,:))), 'o-','Color',[0.10 0.45 0.85],'LineWidth',1.8,'MarkerFaceColor',[0.10 0.45 0.85],'MarkerSize',6,'DisplayName','compensado (self-cal)');
    set(gca,'YScale','log'); xlabel('SNR (dB)'); ylabel('RMSE de DoA (graus)');
    xticks(range_SNR_dB); ylim(ylD); title(methods{m});
    if m==1, legend('Location','best'); end
end
sgtitle(sprintf('DoA por metodo: nao compensado vs compensado  (raio=%.2f\\lambda, %d az x %d real)', ...
    r/lambda, n_angles, n_trials),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'pormetodo_doa_rmse.png'),'Resolution',170);

%% =================== FIGURA 2: ganho de array por metodo ===================
allG = [G_inv(:); G_ideal(:); G_acop(:); G_man(:); G_opt(:)];
ylG = [min(allG)-1, max(allG)+1];
fig = figure('Color','w','Position',[60 60 1200 820]);
for m = 1:nM
    subplot(2,2,m); hold on; grid on;
    plot(range_SNR_dB, G_opt,  '-', 'Color',[0.2 0.2 0.2],'LineWidth',2.6,'DisplayName','otimo (w=Ca)');
    plot(range_SNR_dB, G_man,  'o-','Color',[0.00 0.55 0.15],'LineWidth',2.0,'MarkerFaceColor',[0.00 0.55 0.15],'MarkerSize',6,'DisplayName','manifold (w=b_{hat})');
    plot(range_SNR_dB, G_acop, '--','Color',[0.75 0.10 0.10],'LineWidth',1.6,'DisplayName','acoplado (w=a)');
    plot(range_SNR_dB, G_inv(m,:), 'd-','Color',[0.10 0.25 0.85],'LineWidth',2.0,'MarkerFaceColor',[0.10 0.25 0.85],'MarkerSize',6,'DisplayName','inversao (w=D_m^Ha)');
    set(gca,'YScale','linear'); xlabel('SNR de entrada (dB)'); ylabel('ganho de array (dB)');
    xticks(range_SNR_dB); ylim(ylG); title(methods{m});
    if m==1, legend('Location','best'); end
end
sgtitle(sprintf(['Beamforming por metodo: inverter Chat_m vs manifold  (raio=%.2f\\lambda)\n' ...
    'manifold (independe do metodo) atinge o otimo; inverter fica abaixo'], r/lambda),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'pormetodo_ganho_array.png'),'Resolution',170);

%% ---- Resumo ----
fprintf('\n===== RMSE de DoA compensado (graus), por metodo =====\n');
fprintf('%-8s', 'SNR'); fprintf('%8d', range_SNR_dB); fprintf('\n');
for m=1:nM, fprintf('%-8s', methods{m}); fprintf('%8.3f', squeeze(rmse(m,3,:))); fprintf('\n'); end
fprintf('\n===== Ganho de array inversao (dB), por metodo  (manifold=%.1f dB no fim) =====\n', G_man(end));
for m=1:nM, fprintf('%-8s', methods{m}); fprintf('%8.2f', G_inv(m,:)); fprintf('\n'); end
fprintf('\nFiguras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function c_fin = selfcal_drive(drive_method, X_coupled, b_hat_orig, q, M, r, lambda, ...
                               theta_deg, beta_uca, A_dict, phi_grid_deg, u, maxIter)
% Self-cal dirigido por 'drive_method' (init KW), passo de relaxacao u.
    phi0 = doa_estimate(X_coupled, 'KW', q, r, lambda, beta_uca, A_dict, phi_grid_deg);
    C_hat = estimate_C_circulant_uca(b_hat_orig, utils.steering_vec_uca(M,r,lambda,theta_deg,phi0), M);
    for it = 1:maxIter
        D = safe_inv(C_hat);  Y = D*X_coupled;
        phi = doa_estimate(Y, drive_method, q, r, lambda, beta_uca, A_dict, phi_grid_deg);
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
        case 'DAS'
            R = (X*X')/size(X,2);
            scores = real(sum(conj(A_dict).*(R*A_dict), 1));
            [~,ip] = max(scores); phi_hat = phi_grid_deg(ip);
        case 'CAPON'
            R = (X*X')/size(X,2); R = R + 1e-6*trace(R)/M*eye(M); Rinv = R\eye(M);
            den = real(sum(conj(A_dict).*(Rinv*A_dict),1));
            [~,ip] = max(1./max(den,eps)); phi_hat = phi_grid_deg(ip);
        case 'MUSIC'
            R = (X*X')/size(X,2); R = (R+R')/2;
            [V, Dg] = eig(R); [~, idx] = sort(real(diag(Dg)), 'descend'); V = V(:, idx);
            En = V(:, 2:end); proj = En'*A_dict;
            den = real(sum(conj(proj).*proj, 1));
            [~,ip] = max(1./max(den,eps)); phi_hat = phi_grid_deg(ip);
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
