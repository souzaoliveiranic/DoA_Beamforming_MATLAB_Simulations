% =========================================================================
% DoA_KW_UCA_C_vs_Aperture_N
%
% Estuda o LIMITE RECUPERAVEL da matriz de acoplamento C em funcao de:
%   (1) ABERTURA do array (raio em lambda)  -> condicionamento do problema
%   (2) NUMERO DE AMOSTRAS N (preambulo)     -> ruido de b_hat
%
% Para cada ponto (raio, N) mede-se, em Monte Carlo:
%   - ORACLE 1-DIR (angulo verdadeiro): C por LS circulante de 1 shot. E' o
%     melhor estimador de 1 direcao DADO o angulo, mas NAO e' um piso da
%     metrica de matriz cheia: a estimacao de 1 direcao e' subdeterminada
%     fora do manifold, entao um metodo (angulo estimado) pode bate-lo por
%     acaso na metrica. Serve para ISOLAR o erro de DoA.
%   - ORACLE MULTI-DIR (P_oracle direcoes, angulos verdadeiros): C conjunta
%     bem-determinada -> PISO HONESTO de recuperacao de C (nenhum metodo de
%     1 shot deve bate-lo, a menos de ruido de Monte Carlo).
%   - Self-cal single-shot (4 metodos KW/DAS/Capon/MUSIC, init KW, passo u):
%     o desempenho pratico (angulo estimado).
%
% Metrica de foco: RESIDUO DE COMPENSACAO  ||norm(inv(C_hat) C_true) - I||/sqrt(M)
%   (mede o que importa: quao bem inv(C_hat) desfaz o acoplamento). Frobenius
%   bruto como apoio.
%
% Saidas:
%   (A) vs RAIO  (N fixo): residuo de compensacao e Frobenius, metodos + 2 oracles
%   (B) vs N     (raio fixo): idem
%   (C) GRADE raio x N: oracle multi-dir (piso) e KW (pratico)
% =========================================================================

clear; clc; close all;

% Advertencias ESPERADAS em aberturas pequenas (steering quase paralelos ->
% LS de C mal-condicionado). Sao tratadas (estimativas nao-finitas sao
% DESCARTADAS da media em run_point), entao suprimimos o ruido no console.
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:illConditionedMatrix');
warning('off','estimate_C_circulant_uca:smallAlpha');

%% ---- Parametros fixos ----
M      = 8;  fc = 500e6;  c = 3e8;  lambda = c/fc;  theta_sig_deg = 90;
fs = 288000; Rs = 9600; sps = 30; alpha = 0.3; span = 8; fd = 4.8e3;

maxIter   = 10;
u         = 0.5;                 % passo de relaxacao da atualizacao de C
SNR_fixed = 6;                   % SNR fixa para isolar o efeito de raio e N (dB)
methods   = {'KW','DAS','CAPON','MUSIC'};
nMethods  = numel(methods);
phi_grid_deg = -180:0.5:180;
grid_step = phi_grid_deg(2)-phi_grid_deg(1);
beta_uca  = 2*pi*(0:M-1).'/M;
K         = floor(M/2);

% --- Varreduras ---
range_radius = [0.10 0.15 0.20 0.25];   % em lambda
range_N      = [1050 2100 4200 8400];        % multiplos de sps=30
r_base_idx   = find(abs(range_radius-0.20)<1e-9);  % raio do corte "vs N"
N_base_idx   = find(range_N==2100);                % N do corte "vs raio"
n_real       = 60;               % realizacoes Monte Carlo por ponto
P_oracle     = 5;                % direcoes (ang. verdadeiros) do oracle multi-dir
rng(2026, 'twister');

outDir = fullfile(pwd, 'C vs Aperture-N Graphs');
if ~exist(outDir,'dir'), mkdir(outDir); end

nRad = numel(range_radius); nN = numel(range_N);

%% ---- Pre-computa C_true e dicionario de steering por raio ----
C_true_c   = cell(nRad,1);
normC_c    = zeros(nRad,1);
A_dict_c   = cell(nRad,1);
nGrid = numel(phi_grid_deg);
for ir = 1:nRad
    r = range_radius(ir)*lambda;
    Ct = compute_Ctx_for_R(fc, M, r, 50);
    C_true_c{ir} = Ct;  normC_c(ir) = norm(Ct,'fro');
    A = zeros(M,nGrid);
    for ig = 1:nGrid, A(:,ig) = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_grid_deg(ig)); end
    A_dict_c{ir} = A;
    fprintf('raio=%.2f lambda: ||C_true||=%.3f  cond(C_true)=%.1f\n', range_radius(ir), normC_c(ir), cond(Ct));
end

%% ======================= GRADE raio x N =======================
cres_grid    = zeros(nMethods, nRad, nN);  frob_grid    = zeros(nMethods, nRad, nN);
cres_o1_grid = zeros(nRad, nN);            frob_o1_grid = zeros(nRad, nN);   % oracle 1-dir
cres_om_grid = zeros(nRad, nN);            frob_om_grid = zeros(nRad, nN);   % oracle multi-dir

t0 = tic;
for ir = 1:nRad
    r = range_radius(ir)*lambda;
    for iN = 1:nN
        N = range_N(iN);
        [cm, fm, co1, fo1, com, fom] = run_point(M, r, lambda, theta_sig_deg, SNR_fixed, N, ...
            fs, Rs, sps, alpha, span, fd, phi_grid_deg, grid_step, beta_uca, ...
            A_dict_c{ir}, methods, maxIter, u, C_true_c{ir}, normC_c(ir), n_real, P_oracle);
        cres_grid(:,ir,iN)=cm;  frob_grid(:,ir,iN)=fm;
        cres_o1_grid(ir,iN)=co1; frob_o1_grid(ir,iN)=fo1;
        cres_om_grid(ir,iN)=com; frob_om_grid(ir,iN)=fom;
        fprintf('  raio=%.2f  N=%5d  | oracle1=%.3f  oracleMD=%.3f  KW=%.3f  | %.0fs\n', ...
            range_radius(ir), N, co1, com, cm(1), toc(t0));
    end
end

flr = @(x) max(x, 1e-3);
colorsM = lines(nMethods);
markers_m = {'o-','s-','d-','^-'};

%% ---- (A) vs RAIO (N = range_N(N_base_idx)) ----
iN = N_base_idx;
fig = figure('Color','w','Position',[60 80 1300 520]);
subplot(1,2,1); hold on; grid on;
for im=1:nMethods
    plot(range_radius, flr(squeeze(cres_grid(im,:,iN))), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',methods{im});
end
plot(range_radius, flr(cres_o1_grid(:,iN)), '--p','Color',[.5 .5 .5],'LineWidth',1.6,'MarkerFaceColor',[.5 .5 .5],'MarkerSize',8,'DisplayName','oracle 1-dir (ang. verd.)');
plot(range_radius, flr(cres_om_grid(:,iN)), 'k-p','LineWidth',2.0,'MarkerFaceColor','k','MarkerSize',9,'DisplayName',sprintf('oracle multi-dir P=%d (piso)',P_oracle));
set(gca,'YScale','log'); xlabel('raio (\lambda)'); ylabel('residuo de compensacao'); xticks(range_radius);
title('Residuo de compensacao vs raio'); legend('Location','best');
subplot(1,2,2); hold on; grid on;
for im=1:nMethods
    plot(range_radius, flr(squeeze(frob_grid(im,:,iN))), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',methods{im});
end
plot(range_radius, flr(frob_o1_grid(:,iN)), '--p','Color',[.5 .5 .5],'LineWidth',1.6,'MarkerFaceColor',[.5 .5 .5],'MarkerSize',8,'DisplayName','oracle 1-dir');
plot(range_radius, flr(frob_om_grid(:,iN)), 'k-p','LineWidth',2.0,'MarkerFaceColor','k','MarkerSize',9,'DisplayName',sprintf('oracle multi-dir P=%d',P_oracle));
set(gca,'YScale','log'); xlabel('raio (\lambda)'); ylabel('Frobenius bruto'); xticks(range_radius);
title('Frobenius vs raio'); legend('Location','best');
sgtitle(sprintf('Influencia da ABERTURA  (N=%d, SNR=%+d dB, %d realizacoes, u=%.2f)', ...
    range_N(iN), SNR_fixed, n_real, u),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'C_vs_raio.png'),'Resolution',170);

%% ---- (B) vs N (raio = range_radius(r_base_idx)) ----
ir = r_base_idx;
fig = figure('Color','w','Position',[60 80 1300 520]);
subplot(1,2,1); hold on; grid on;
for im=1:nMethods
    plot(range_N, flr(squeeze(cres_grid(im,ir,:))), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',methods{im});
end
plot(range_N, flr(cres_o1_grid(ir,:)), '--p','Color',[.5 .5 .5],'LineWidth',1.6,'MarkerFaceColor',[.5 .5 .5],'MarkerSize',8,'DisplayName','oracle 1-dir (ang. verd.)');
plot(range_N, flr(cres_om_grid(ir,:)), 'k-p','LineWidth',2.0,'MarkerFaceColor','k','MarkerSize',9,'DisplayName',sprintf('oracle multi-dir P=%d (piso)',P_oracle));
set(gca,'YScale','log','XScale','log'); xlabel('N (amostras)'); ylabel('residuo de compensacao'); xticks(range_N);
title('Residuo de compensacao vs N'); legend('Location','best');
subplot(1,2,2); hold on; grid on;
for im=1:nMethods
    plot(range_N, flr(squeeze(frob_grid(im,ir,:))), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',methods{im});
end
plot(range_N, flr(frob_o1_grid(ir,:)), '--p','Color',[.5 .5 .5],'LineWidth',1.6,'MarkerFaceColor',[.5 .5 .5],'MarkerSize',8,'DisplayName','oracle 1-dir');
plot(range_N, flr(frob_om_grid(ir,:)), 'k-p','LineWidth',2.0,'MarkerFaceColor','k','MarkerSize',9,'DisplayName',sprintf('oracle multi-dir P=%d',P_oracle));
set(gca,'YScale','log','XScale','log'); xlabel('N (amostras)'); ylabel('Frobenius bruto'); xticks(range_N);
title('Frobenius vs N'); legend('Location','best');
sgtitle(sprintf('Influencia de N  (raio=%.2f\\lambda, SNR=%+d dB, %d realizacoes, u=%.2f)', ...
    range_radius(ir), SNR_fixed, n_real, u),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'C_vs_N.png'),'Resolution',170);

%% ---- (C) GRADE raio x N: piso (oracle multi-dir) e pratico (KW) ----
rad_colors = cool(nRad);
fig = figure('Color','w','Position',[60 80 1300 520]);
subplot(1,2,1); hold on; grid on;
for ir=1:nRad
    plot(range_N, flr(cres_om_grid(ir,:)), 'o-','Color',rad_colors(ir,:),'LineWidth',1.8, ...
        'MarkerFaceColor',rad_colors(ir,:),'MarkerSize',8,'DisplayName',sprintf('%.2f\\lambda',range_radius(ir)));
end
set(gca,'YScale','log','XScale','log'); xlabel('N (amostras)'); ylabel('residuo de compensacao'); xticks(range_N);
title(sprintf('PISO: oracle multi-dir P=%d, compRes vs N por raio',P_oracle)); legend('Location','best');
subplot(1,2,2); hold on; grid on;
iKW = find(strcmp(methods,'KW'));
for ir=1:nRad
    plot(range_N, flr(squeeze(cres_grid(iKW,ir,:))), 'o-','Color',rad_colors(ir,:),'LineWidth',1.8, ...
        'MarkerFaceColor',rad_colors(ir,:),'MarkerSize',8,'DisplayName',sprintf('%.2f\\lambda',range_radius(ir)));
end
set(gca,'YScale','log','XScale','log'); xlabel('N (amostras)'); ylabel('residuo de compensacao'); xticks(range_N);
title('Self-cal KW (pratico): compRes vs N por raio'); legend('Location','best');
sgtitle(sprintf('Limite recuperavel de C: raio x N  (SNR=%+d dB, %d realizacoes)', SNR_fixed, n_real),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'C_grade_raio_x_N.png'),'Resolution',170);

%% ---- Resumo numerico ----
fprintf('\n===== RESIDUO DE COMPENSACAO (oracle multi-dir P=%d) — grade raio x N =====\n', P_oracle);
fprintf('%-8s', 'raio\\N'); fprintf('%9d', range_N); fprintf('\n');
for ir=1:nRad
    fprintf('%-8.2f', range_radius(ir)); fprintf('%9.3f', cres_om_grid(ir,:)); fprintf('\n');
end
fprintf('\nFiguras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function [cres_m, frob_m, cres_o1, frob_o1, cres_om, frob_om] = run_point(M, r, lambda, ...
    theta, SNR, N, fs, Rs, sps, alpha, span, fd, phi_grid, grid_step, beta, A_dict, ...
    methods, maxIter, u, C_true, normCtrue, n_real, P_oracle)
    nM = numel(methods);
    % Somas e CONTAGENS de validos (estimativas nao-finitas sao descartadas).
    sc_m=zeros(nM,1); sf_m=zeros(nM,1); cnt_m=zeros(nM,1);
    sc_o1=0; sf_o1=0; cnt_o1=0;
    sc_om=0; sf_om=0; cnt_om=0;
    Bbuf = zeros(M, P_oracle); Abuf = zeros(M, P_oracle); fillb = 0;
    nGroups = floor(n_real/P_oracle);
    for run = 1:n_real
        phi_true = round((-180 + 360*rand)/grid_step)*grid_step;
        [~, q, ~, Xsig, ~, Xn] = utils.simulate_fsk_data_uca(M, r, lambda, ...
            phi_true, phi_true+90, theta, theta, SNR, -100, N, fs, Rs, sps, alpha, span, fd);
        q = q(:);  Xc = C_true*Xsig + Xn;  bhat = Xc*conj(q)/(q'*q);
        a_true = utils.steering_vec_uca(M, r, lambda, theta, phi_true);

        % --- oracle 1-dir ---
        Co1 = estimate_C_circulant_uca(bhat, a_true, M);
        [sc_o1, sf_o1, cnt_o1] = accum(sc_o1, sf_o1, cnt_o1, Co1, C_true, normCtrue);

        % --- buffer para oracle multi-dir (reaproveita shots; sem custo extra) ---
        fillb = fillb + 1;  Bbuf(:,fillb) = bhat;  Abuf(:,fillb) = a_true;
        if fillb == P_oracle
            Com = estimate_C_multidir(Bbuf, Abuf, M, 30, 1e-10);
            [sc_om, sf_om, cnt_om] = accum(sc_om, sf_om, cnt_om, Com, C_true, normCtrue);
            fillb = 0;
        end

        % --- self-cal por metodo (1 shot) ---
        for im = 1:nM
            [frob_it, cres_it] = selfcal_one_damped(methods{im}, Xc, bhat, q, ...
                M, r, lambda, theta, beta, A_dict, phi_grid, C_true, normCtrue, u, maxIter);
            cr = cres_it(end); fr = frob_it(end);
            if isfinite(cr) && isfinite(fr)
                sc_m(im)=sc_m(im)+cr; sf_m(im)=sf_m(im)+fr; cnt_m(im)=cnt_m(im)+1;
            end
        end
    end
    cres_m = sc_m./max(cnt_m,1);  frob_m = sf_m./max(cnt_m,1);
    cres_o1= sc_o1/max(cnt_o1,1); frob_o1= sf_o1/max(cnt_o1,1);
    cres_om= sc_om/max(cnt_om,1); frob_om= sf_om/max(cnt_om,1);

    % Aviso de descarte (estimativas degeneradas em abertura pequena)
    ndrop_o1 = n_real - cnt_o1;  ndrop_om = nGroups - cnt_om;  ndrop_m = max(n_real - cnt_m);
    if ndrop_o1>0 || ndrop_om>0 || ndrop_m>0
        fprintf('     [descartes] oracle1=%d/%d  oracleMD=%d/%d  metodo(pior)=%d/%d\n', ...
            ndrop_o1, n_real, ndrop_om, nGroups, ndrop_m, n_real);
    end
end

function [sc, sf, cnt] = accum(sc, sf, cnt, C_hat, C_true, normCtrue)
% Acumula residuo de compensacao e Frobenius SE a estimativa for finita.
    cr = comp_residual(C_hat, C_true);
    fr = norm(C_hat - C_true,'fro')/normCtrue;
    if isfinite(cr) && isfinite(fr)
        sc = sc + cr;  sf = sf + fr;  cnt = cnt + 1;
    end
end

function [frob_it, cres_it] = selfcal_one_damped(method, X_coupled, b_hat_orig, q, ...
    M, r, lambda, theta_deg, beta_uca, A_dict, phi_grid_deg, C_true, normCtrue, u, maxIter)
    phi0  = doa_estimate(X_coupled, 'KW', q, r, lambda, beta_uca, A_dict, phi_grid_deg);
    a0    = utils.steering_vec_uca(M, r, lambda, theta_deg, phi0);
    C_hat = estimate_C_circulant_uca(b_hat_orig, a0, M);
    frob_it=zeros(maxIter,1); cres_it=zeros(maxIter,1);
    for it = 1:maxIter
        D   = safe_inv(C_hat);
        Y   = D * X_coupled;
        phi = doa_estimate(Y, method, q, r, lambda, beta_uca, A_dict, phi_grid_deg);
        a_hat = utils.steering_vec_uca(M, r, lambda, theta_deg, phi);
        C_ls  = estimate_C_circulant_uca(b_hat_orig, a_hat, M);
        C_hat = (1-u)*C_hat + u*C_ls;
        frob_it(it) = norm(C_hat - C_true,'fro')/normCtrue;
        cres_it(it) = comp_residual(C_hat, C_true);
    end
end

function rr = comp_residual(C_hat, C_true)
% Residuo de compensacao efetiva. Retorna NaN (estimativa INVALIDA) se C_hat
% nao for finita -> sera DESCARTADA da media (em vez de poluir com um valor).
    if ~all(isfinite(C_hat(:))), rr = NaN; return; end
    M = size(C_hat,1);
    E = safe_inv(C_hat) * C_true;
    g = trace(E)/M; if abs(g) > eps, E = E / g; end
    if all(isfinite(E(:))), rr = norm(E - eye(M),'fro')/sqrt(M);
    else,                   rr = NaN; end
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
