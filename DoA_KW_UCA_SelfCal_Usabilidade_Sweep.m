% =========================================================================
% DoA_KW_UCA_SelfCal_Usabilidade_Sweep
%
% Varredura de ABERTURA (raio) e de N para responder: "em que regime, se
% algum, vale a pena ESTIMAR C e INVERTER (desacoplar)?" vs nao compensar vs
% calibrar o manifold (w=b_hat, sem inverter).
%
% Metrica: perda de SNR de array (dB) = 10log10(G_id/G), MENOR = melhor.
%   G_dec = |a^H D C a|^2/(a^H D D^H a),  D=inv(Chat), w=a   (desacopla)
%   G_nc  = |a^H C a|^2/(a^H a)                              (sem comp.)
%   G_man = |b_hat^H C a|^2/(b_hat^H b_hat)                  (manifold, SEM inverter)
%   G_opt = ||C a||^2                                        (matched ao canal)
%
% (1) perda vs RAIO  (N e SNR fixos)
% (2) perda vs N     (raio e SNR fixos)
% A decisao le-se direto: onde a curva de desacoplamento (oracle = melhor
% estimacao possivel) fica ABAIXO de sem-comp E de manifold.
% =========================================================================

clear; clc; close all;

warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:illConditionedMatrix');
warning('off','estimate_C_circulant_uca:smallAlpha');

%% ---- Parametros fixos ----
M = 8; fc = 500e6; c = 3e8; lambda = c/fc; theta_sig_deg = 90;
fs = 288000; Rs = 9600; sps = 30; alpha = 0.3; span = 8; fd = 4.8e3;

maxIter   = 10;
u         = 0.5;
SNR_fixed = 6;                  % SNR fixa para isolar raio/N (dB)
methods   = {'KW','DAS','CAPON','MUSIC'};
nMethods  = numel(methods);
phi_grid_deg = -180:0.5:180;
grid_step = phi_grid_deg(2)-phi_grid_deg(1);
beta_uca  = 2*pi*(0:M-1).'/M;

% --- Varreduras ---
range_radius = [0.10 0.15 0.20 0.25 0.30 0.40 0.50];   % em lambda
range_N      = [1050 2100 4200 8400];
r_base       = 0.20;            % raio do corte "vs N"
N_base       = 8400;            % N do corte "vs raio"
n_real       = 100;             % realizacoes Monte Carlo por ponto
rng(2026, 'twister');

outDir = fullfile(pwd, 'Usabilidade Graphs');
if ~exist(outDir,'dir'), mkdir(outDir); end

nRad = numel(range_radius); nN = numel(range_N);

%% ======================= (1) VARREDURA DE RAIO =======================
fprintf('=== Varredura de RAIO (N=%d, SNR=%+d dB) ===\n', N_base, SNR_fixed);
Gid_r=zeros(1,nRad); Gnc_r=zeros(1,nRad); Gman_r=zeros(1,nRad); Gopt_r=zeros(1,nRad);
Gdec_r=zeros(nMethods,nRad); Gdec_o_r=zeros(1,nRad);
t0 = tic;
for ir = 1:nRad
    r = range_radius(ir)*lambda;
    C_true = compute_Ctx_for_R(fc, M, r, 50);  normCtrue = norm(C_true,'fro');
    A_dict = build_dict(M, r, lambda, theta_sig_deg, phi_grid_deg);
    R = run_point(M, r, lambda, theta_sig_deg, SNR_fixed, N_base, fs,Rs,sps,alpha,span,fd, ...
        phi_grid_deg, grid_step, beta_uca, A_dict, methods, maxIter, u, C_true, normCtrue, n_real);
    Gid_r(ir)=R.Gid; Gnc_r(ir)=R.Gnc; Gman_r(ir)=R.Gman; Gopt_r(ir)=R.Gopt;
    Gdec_r(:,ir)=R.Gdec; Gdec_o_r(ir)=R.Gdec_o;
    fprintf('  raio=%.2f  cond=%.0f  | %.0fs\n', range_radius(ir), cond(C_true), toc(t0));
end

%% ======================= (2) VARREDURA DE N =======================
fprintf('\n=== Varredura de N (raio=%.2f lambda, SNR=%+d dB) ===\n', r_base, SNR_fixed);
rN = r_base*lambda;  C_trueN = compute_Ctx_for_R(fc, M, rN, 50);  normCtrueN = norm(C_trueN,'fro');
A_dictN = build_dict(M, rN, lambda, theta_sig_deg, phi_grid_deg);
Gid_N=zeros(1,nN); Gnc_N=zeros(1,nN); Gman_N=zeros(1,nN); Gopt_N=zeros(1,nN);
Gdec_N=zeros(nMethods,nN); Gdec_o_N=zeros(1,nN);
for iN = 1:nN
    Nv = range_N(iN);
    R = run_point(M, rN, lambda, theta_sig_deg, SNR_fixed, Nv, fs,Rs,sps,alpha,span,fd, ...
        phi_grid_deg, grid_step, beta_uca, A_dictN, methods, maxIter, u, C_trueN, normCtrueN, n_real);
    Gid_N(iN)=R.Gid; Gnc_N(iN)=R.Gnc; Gman_N(iN)=R.Gman; Gopt_N(iN)=R.Gopt;
    Gdec_N(:,iN)=R.Gdec; Gdec_o_N(iN)=R.Gdec_o;
    fprintf('  N=%5d  | %.0fs\n', Nv, toc(t0));
end

%% ---- Perda em dB ----
LdB = @(Gid,G) 10*log10(Gid ./ G);
% raio
loss_nc_r=LdB(Gid_r,Gnc_r); loss_man_r=LdB(Gid_r,Gman_r); loss_opt_r=LdB(Gid_r,Gopt_r);
loss_dec_r=10*log10(Gid_r./Gdec_r); loss_dec_o_r=LdB(Gid_r,Gdec_o_r);
% N
loss_nc_N=LdB(Gid_N,Gnc_N); loss_man_N=LdB(Gid_N,Gman_N); loss_opt_N=LdB(Gid_N,Gopt_N);
loss_dec_N=10*log10(Gid_N./Gdec_N); loss_dec_o_N=LdB(Gid_N,Gdec_o_N);

colorsM = lines(nMethods); markers_m = {'o-','s-','d-','^-'};

%% ---- Figura (1): perda vs RAIO ----
fig = figure('Color','w','Position',[60 80 1000 620]); hold on; grid on;
for im=1:nMethods
    plot(range_radius, loss_dec_r(im,:), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',sprintf('desacopla %s',methods{im}));
end
plot(range_radius, loss_dec_o_r,'k--p','LineWidth',1.8,'MarkerFaceColor','k','MarkerSize',8,'DisplayName','desacopla ORACLE');
plot(range_radius, loss_nc_r, ':','Color',[.45 .45 .45],'LineWidth',2.4,'DisplayName','SEM compensacao');
plot(range_radius, loss_man_r,'-.','Color',[0 0.55 0],'LineWidth',2.4,'DisplayName','MANIFOLD w=b_{hat}');
plot(range_radius, loss_opt_r,'-','Color',[0.3 0.3 0.3],'LineWidth',1.0,'DisplayName','otimo (matched)');
yline(0,'r--','LineWidth',1.0,'DisplayName','ideal (sem acoplamento)');
set(gca,'YDir','reverse'); xlabel('raio (\lambda)'); ylabel('perda de SNR de array (dB) [menor=melhor]');
xticks(range_radius); title(sprintf('Usabilidade vs ABERTURA  (N=%d, SNR=%+d dB, %d realizacoes, u=%.2f)', N_base, SNR_fixed, n_real, u));
legend('Location','eastoutside');
exportgraphics(fig, fullfile(outDir,'usabilidade_vs_raio.png'),'Resolution',180);
matlab2tikz(fullfile(outDir,'usabilidade_vs_raio.tex'), 'width','\figurewidth','height','\figureheight');

%% ---- Figura (2): perda vs N ----
fig = figure('Color','w','Position',[60 80 1000 620]); hold on; grid on;
for im=1:nMethods
    plot(range_N, loss_dec_N(im,:), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',sprintf('desacopla %s',methods{im}));
end
plot(range_N, loss_dec_o_N,'k--p','LineWidth',1.8,'MarkerFaceColor','k','MarkerSize',8,'DisplayName','desacopla ORACLE');
plot(range_N, loss_nc_N, ':','Color',[.45 .45 .45],'LineWidth',2.4,'DisplayName','SEM compensacao');
plot(range_N, loss_man_N,'-.','Color',[0 0.55 0],'LineWidth',2.4,'DisplayName','MANIFOLD w=b_{hat}');
plot(range_N, loss_opt_N,'-','Color',[0.3 0.3 0.3],'LineWidth',1.0,'DisplayName','otimo (matched)');
yline(0,'r--','LineWidth',1.0,'DisplayName','ideal (sem acoplamento)');
set(gca,'YDir','reverse','XScale','log'); xlabel('N (amostras)'); ylabel('perda de SNR de array (dB) [menor=melhor]');
xticks(range_N); title(sprintf('Usabilidade vs N  (raio=%.2f\\lambda, SNR=%+d dB, %d realizacoes, u=%.2f)', r_base, SNR_fixed, n_real, u));
legend('Location','eastoutside');
exportgraphics(fig, fullfile(outDir,'usabilidade_vs_N.png'),'Resolution',180);
matlab2tikz(fullfile(outDir,'usabilidade_vs_N.tex'), 'width','\figurewidth','height','\figureheight');

%% ---- Resumo numerico (raio) ----
fprintf('\n===== Perda de SNR (dB) vs RAIO (N=%d, SNR=%+d dB) =====\n', N_base, SNR_fixed);
fprintf('%-8s %9s %9s %9s %9s\n','raio','semComp','manifold','dec-KW','dec-ORC');
for ir=1:nRad
    fprintf('%-8.2f %+9.2f %+9.2f %+9.2f %+9.2f\n', range_radius(ir), ...
        loss_nc_r(ir), loss_man_r(ir), loss_dec_r(1,ir), loss_dec_o_r(ir));
end
fprintf('\n(desacopla < semComp E < manifold -> vale inverter C; senao, nao vale)\n');
fprintf('Figuras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function A = build_dict(M, r, lambda, theta, phi_grid)
    A = zeros(M, numel(phi_grid));
    for ig = 1:numel(phi_grid), A(:,ig) = utils.steering_vec_uca(M, r, lambda, theta, phi_grid(ig)); end
end

function R = run_point(M, r, lambda, theta, SNR, N, fs,Rs,sps,alpha,span,fd, ...
    phi_grid, grid_step, beta, A_dict, methods, maxIter, u, C_true, normCtrue, n_real)
    nM = numel(methods);
    sGid=0; sGnc=0; sGman=0; sGopt=0; cd=0;
    sGdec=zeros(nM,1); cm=zeros(nM,1);
    sGdec_o=0; co=0;
    for run = 1:n_real
        phi = round((-180 + 360*rand)/grid_step)*grid_step;
        a   = utils.steering_vec_uca(M, r, lambda, theta, phi);  Ca = C_true*a;
        [~, q, ~, Xsig, ~, Xn] = utils.simulate_fsk_data_uca(M, r, lambda, ...
            phi, phi+90, theta, theta, SNR, -100, N, fs, Rs, sps, alpha, span, fd);
        q = q(:);  Xc = C_true*Xsig + Xn;  bhat = Xc*conj(q)/(q'*q);

        sGid  = sGid  + real(a'*a);
        sGnc  = sGnc  + abs(a'*Ca)^2/real(a'*a);
        sGopt = sGopt + real(Ca'*Ca);
        sGman = sGman + abs(bhat'*Ca)^2/max(real(bhat'*bhat),eps);
        cd = cd + 1;

        [~,~,~,corc] = selfcal_one_damped('ORACLE', Xc, bhat, q, phi, M, r, lambda, theta, ...
            beta, A_dict, phi_grid, C_true, normCtrue, u, maxIter);
        Gd = usability(reconstruct_C(corc,M), C_true, a);
        if isfinite(Gd), sGdec_o = sGdec_o + Gd; co = co + 1; end

        for im = 1:nM
            [~,~,~,cf] = selfcal_one_damped(methods{im}, Xc, bhat, q, phi, M, r, lambda, theta, ...
                beta, A_dict, phi_grid, C_true, normCtrue, u, maxIter);
            Gd = usability(reconstruct_C(cf,M), C_true, a);
            if isfinite(Gd), sGdec(im) = sGdec(im) + Gd; cm(im) = cm(im) + 1; end
        end
    end
    R.Gid=sGid/cd; R.Gnc=sGnc/cd; R.Gman=sGman/cd; R.Gopt=sGopt/cd;
    R.Gdec=sGdec./max(cm,1); R.Gdec_o=sGdec_o/max(co,1);
end

function Gdec = usability(C_hat, C_true, a)
    if ~all(isfinite(C_hat(:))), Gdec=NaN; return; end
    D   = safe_inv(C_hat);
    DCa = D*(C_true*a);
    sig = a'*DCa;
    noi = real(a'*(D*D')*a);
    if noi>eps && isfinite(noi), Gdec = abs(sig)^2/noi; else, Gdec = NaN; end
end

function [frob_it, cres_it, derr_it, c_fin] = selfcal_one_damped(method, X_coupled, ...
    b_hat_orig, q, phi_true, M, r, lambda, theta_deg, beta_uca, A_dict, phi_grid_deg, ...
    C_true, normCtrue, u, maxIter)
    phi0  = doa_estimate(X_coupled, 'KW', q, r, lambda, beta_uca, A_dict, phi_grid_deg);
    a0    = utils.steering_vec_uca(M, r, lambda, theta_deg, phi0);
    C_hat = estimate_C_circulant_uca(b_hat_orig, a0, M);
    frob_it=zeros(maxIter,1); cres_it=zeros(maxIter,1); derr_it=zeros(maxIter,1);
    for it = 1:maxIter
        D   = safe_inv(C_hat);
        Y   = D * X_coupled;
        if strcmpi(method,'ORACLE'), phi = phi_true;
        else, phi = doa_estimate(Y, method, q, r, lambda, beta_uca, A_dict, phi_grid_deg); end
        a_hat = utils.steering_vec_uca(M, r, lambda, theta_deg, phi);
        C_ls  = estimate_C_circulant_uca(b_hat_orig, a_hat, M);
        C_hat = (1-u)*C_hat + u*C_ls;
        frob_it(it) = norm(C_hat - C_true,'fro')/normCtrue;
        cres_it(it) = norm(safe_inv(C_hat)*C_true - eye(M),'fro')/sqrt(M);
        derr_it(it) = abs(wrapTo180(phi - phi_true));
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
