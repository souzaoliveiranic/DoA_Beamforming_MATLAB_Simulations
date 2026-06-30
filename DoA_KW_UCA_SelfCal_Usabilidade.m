% =========================================================================
% DoA_KW_UCA_SelfCal_Usabilidade
%
% Derivado de DoA_KW_UCA_SelfCal_Iterativo_Completo. Responde: "a matriz de
% acoplamento estimada tem USABILIDADE como matriz de DESACOPLAMENTO?" — ou
% seja, adianta estimar C para usar D=inv(C) no processamento futuro?
%
% Metrica decisiva: GANHO DE ARRAY (SNR de saida) apos o processamento, que
% funde os DOIS modos de falha: acoplamento residual E amplificacao de ruido.
% Com 'a' normalizado (||a||=1) e peso w=a (aponta na direcao verdadeira):
%
%   (ideal, sem acoplamento)         G_id  = |a^H a|^2 / (a^H a) = 1
%   (sem comp., w=a no dado acopl.)  G_nc  = |a^H C a|^2 / (a^H a)
%   (DESACOPLA D=inv(Chat), w=a)     G_dec = |a^H D C a|^2 / (a^H D D^H a)
%   (MANIFOLD calibrado, w=b_hat)    G_man = |b_hat^H C a|^2 / (b_hat^H b_hat)
%   (otimo: matched ao canal real)   G_opt = ||C a||^2
%
% Perda de SNR (dB) = 10log10(G_id / G_X).  Menor = melhor.  0 dB = ideal.
%
% (a) Perda de SNR vs SNR: ideal / sem-comp / desacoplado (4 metodos) / oracle.
% (b) Decomposicao: amplificacao de ruido trace(DD^H)/M e fidelidade de steering.
% (c) Manifold calibrado (w=b_hat, SEM inverter / SEM estimar C) + teto G_opt.
%
% IMPORTANTE: o manifold (c) usa o b_hat MEDIDO como peso -> nao precisa
% estimar C. Se ele vencer o desacoplamento, a conclusao e' que estimar C
% para inverter NAO adianta nesse cenario.
% =========================================================================

clear; clc; close all;

warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:illConditionedMatrix');
warning('off','estimate_C_circulant_uca:smallAlpha');

%% ---- Parametros do array ----
M      = 8;  fc = 500e6;  c = 3e8;  lambda = c/fc;
r      = 0.15*lambda;       % <-- raio (escolhivel)
theta_sig_deg = 90;

%% ---- Parametros do sinal ----
fs = 288000; N = 3000; Rs = 9600; sps = 30; alpha = 0.3; span = 8; fd = 4.8e3;

%% ---- Parametros do experimento ----
maxIter      = 6;
u            = 0.5;                 % <-- passo de relaxacao ESCOLHIVEL
range_SNR_dB = -12:3:12; %[-6, 0, 6, 12];
methods      = {'KW','DAS','CAPON','MUSIC'};
nMethods     = numel(methods);
phi_grid_deg = -180:0.5:180;
grid_step    = phi_grid_deg(2) - phi_grid_deg(1);
beta_uca     = 2*pi*(0:M-1).'/M;
K            = floor(M/2);

% --- Monte Carlo ---
n_angles = 200;  n_trials = 5;  n_real = n_angles*n_trials;
rng(2026, 'twister');
phi_set  = round( (-180 + 360*rand(1, n_angles)) / grid_step ) * grid_step;

outDir = fullfile(pwd, 'Usabilidade Graphs');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---- C_true e dicionario ----
C_true    = compute_Ctx_for_R(fc, M, r, 50);
normCtrue = norm(C_true, 'fro');
fprintf('C_true (raio=%.2f lambda, ||C_true||_F=%.4f, cond=%.1f), u=%.2f\n', ...
    r/lambda, normCtrue, cond(C_true), u);
nGrid = numel(phi_grid_deg);
A_dict = zeros(M, nGrid);
for ig = 1:nGrid, A_dict(:,ig) = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_grid_deg(ig)); end

%% ---- Acumuladores (somas de ganho LINEAR + contagem de validos) ----
nSNR = numel(range_SNR_dB);
% por SNR (independem do metodo)
s_Gid=zeros(1,nSNR); s_Gnc=zeros(1,nSNR); s_Gopt=zeros(1,nSNR); s_Gman=zeros(1,nSNR); cnt_dir=zeros(1,nSNR);
% por metodo x SNR (desacoplamento)
s_Gdec=zeros(nMethods,nSNR); s_namp=zeros(nMethods,nSNR); s_fid=zeros(nMethods,nSNR); cnt_m=zeros(nMethods,nSNR);
% oracle x SNR
s_Gdec_o=zeros(1,nSNR); s_namp_o=zeros(1,nSNR); s_fid_o=zeros(1,nSNR); cnt_o=zeros(1,nSNR);

t0 = tic;
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fprintf('\n===== SNR = %+d dB  (u=%.2f) =====\n', SNR_dB, u);
    for ia = 1:n_angles
        phi_true = phi_set(ia);
        a = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_true);   % ||a||=1
        Ca = C_true*a;
        for itr = 1:n_trials
            [~, q, ~, Xsig, ~, Xn] = utils.simulate_fsk_data_uca(M, r, lambda, ...
                phi_true, phi_true+90, theta_sig_deg, theta_sig_deg, ...
                SNR_dB, -100, N, fs, Rs, sps, alpha, span, fd);
            q = q(:);
            X_coupled  = C_true*Xsig + Xn;
            b_hat      = X_coupled*conj(q)/(q'*q);

            % --- ganhos que NAO dependem de Chat ---
            s_Gid(iSNR)  = s_Gid(iSNR)  + real(a'*a);
            s_Gnc(iSNR)  = s_Gnc(iSNR)  + abs(a'*Ca)^2/real(a'*a);
            s_Gopt(iSNR) = s_Gopt(iSNR) + real(Ca'*Ca);
            s_Gman(iSNR) = s_Gman(iSNR) + abs(b_hat'*Ca)^2/max(real(b_hat'*b_hat),eps);
            cnt_dir(iSNR)= cnt_dir(iSNR)+ 1;

            % --- ORACLE (C com angulo verdadeiro, mesmo laco/u) ---
            [~,~,~,c_orc] = selfcal_one_damped('ORACLE', X_coupled, b_hat, q, phi_true, ...
                M, r, lambda, theta_sig_deg, beta_uca, A_dict, phi_grid_deg, C_true, normCtrue, u, maxIter);
            [Gd,na,fi] = usability(reconstruct_C(c_orc,M), C_true, a);
            if isfinite(Gd)&&isfinite(na)&&isfinite(fi)
                s_Gdec_o(iSNR)=s_Gdec_o(iSNR)+Gd; s_namp_o(iSNR)=s_namp_o(iSNR)+na; s_fid_o(iSNR)=s_fid_o(iSNR)+fi; cnt_o(iSNR)=cnt_o(iSNR)+1;
            end

            % --- metodos (Chat estimado) ---
            for im = 1:nMethods
                [~,~,~,c_fin] = selfcal_one_damped(methods{im}, X_coupled, b_hat, q, phi_true, ...
                    M, r, lambda, theta_sig_deg, beta_uca, A_dict, phi_grid_deg, C_true, normCtrue, u, maxIter);
                [Gd,na,fi] = usability(reconstruct_C(c_fin,M), C_true, a);
                if isfinite(Gd)&&isfinite(na)&&isfinite(fi)
                    s_Gdec(im,iSNR)=s_Gdec(im,iSNR)+Gd; s_namp(im,iSNR)=s_namp(im,iSNR)+na; s_fid(im,iSNR)=s_fid(im,iSNR)+fi; cnt_m(im,iSNR)=cnt_m(im,iSNR)+1;
                end
            end
        end
        fprintf('  phi=%+6.1f deg | %.0fs\n', phi_true, toc(t0));
    end
end

%% ---- Consolidacao (media linear -> perda em dB) ----
mGid = s_Gid./cnt_dir;  mGnc = s_Gnc./cnt_dir;  mGopt = s_Gopt./cnt_dir;  mGman = s_Gman./cnt_dir;
mGdec = s_Gdec./max(cnt_m,1);  mnamp = s_namp./max(cnt_m,1);  mfid = s_fid./max(cnt_m,1);
mGdec_o = s_Gdec_o./max(cnt_o,1);  mnamp_o = s_namp_o./max(cnt_o,1);  mfid_o = s_fid_o./max(cnt_o,1);

lossdB     = @(G) 10*log10(mGid ./ G);     % perda vs ideal (linha por SNR)
loss_nc    = lossdB(mGnc);
loss_man   = lossdB(mGman);
loss_opt   = lossdB(mGopt);
loss_dec   = 10*log10(mGid ./ mGdec);      % nMethods x nSNR (broadcast mGid 1xnSNR)
loss_dec_o = lossdB(mGdec_o);

colorsM = lines(nMethods);
markers_m = {'o-','s-','d-','^-'};

%% ---- (a)+(c) Perda de SNR de array vs SNR ----
fig = figure('Color','w','Position',[80 80 980 600]); hold on; grid on;
% desacoplamento por metodo
for im = 1:nMethods
    plot(range_SNR_dB, loss_dec(im,:), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8, ...
        'DisplayName',sprintf('desacopla %s',methods{im}));
end
plot(range_SNR_dB, loss_dec_o, 'k--p','LineWidth',1.8,'MarkerFaceColor','k','MarkerSize',8,'DisplayName','desacopla ORACLE');
plot(range_SNR_dB, loss_nc,  ':','Color',[.45 .45 .45],'LineWidth',2.2,'DisplayName','SEM compensacao (w=a)');
plot(range_SNR_dB, loss_man, '-.','Color',[0 0.55 0],'LineWidth',2.2,'DisplayName','MANIFOLD w=b_{hat} (sem inverter)');
plot(range_SNR_dB, loss_opt, 'r-','LineWidth',1.2,'DisplayName','otimo (matched ao canal)');
yline(0,'r--','LineWidth',1.0,'DisplayName','ideal (sem acoplamento)');
set(gca,'YDir','reverse');   % menor perda em cima
xlabel('SNR (dB)'); ylabel('perda de SNR de array (dB)  [menor = melhor]');
xticks(range_SNR_dB); xlim([min(range_SNR_dB)-1, max(range_SNR_dB)+1]);
title(sprintf(['Usabilidade do desacoplamento: perda de SNR vs SNR\n' ...
    'u=%.2f, raio=%.2f\\lambda, %d realizacoes'], u, r/lambda, n_real));
legend('Location','eastoutside');
exportgraphics(fig, fullfile(outDir,'usabilidade_perda_SNR_vs_snr.png'),'Resolution',180);
matlab2tikz(fullfile(outDir,'usabilidade_perda_SNR_vs_snr.tex'), 'width','\figurewidth','height','\figureheight');

%% ---- (b) Decomposicao: amplificacao de ruido e fidelidade de steering ----
fig = figure('Color','w','Position',[60 80 1300 500]);
subplot(1,2,1); hold on; grid on;
for im = 1:nMethods
    plot(range_SNR_dB, 10*log10(mnamp(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',methods{im});
end
plot(range_SNR_dB, 10*log10(mnamp_o), 'k--p','LineWidth',1.8,'MarkerFaceColor','k','MarkerSize',8,'DisplayName','oracle');
xlabel('SNR (dB)'); ylabel('amplificacao de ruido  10log_{10}(trace(DD^H)/M)  (dB)');
xticks(range_SNR_dB); title('Amplificacao de ruido do desacoplador  D=inv(C_{hat})'); legend('Location','best');

subplot(1,2,2); hold on; grid on;
for im = 1:nMethods
    plot(range_SNR_dB, mfid(im,:), markers_m{im}, 'Color',colorsM(im,:), ...
        'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',methods{im});
end
plot(range_SNR_dB, mfid_o, 'k--p','LineWidth',1.8,'MarkerFaceColor','k','MarkerSize',8,'DisplayName','oracle');
xlabel('SNR (dB)'); ylabel('fidelidade de steering  |a^H DCa|^2/(||a||^2||DCa||^2)');
xticks(range_SNR_dB); ylim([0 1.05]); title('Fidelidade de steering (1 = sem acoplamento residual)'); legend('Location','best');
sgtitle(sprintf('Decomposicao da usabilidade  (u=%.2f, raio=%.2f\\lambda)', u, r/lambda),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'usabilidade_decomposicao.png'),'Resolution',170);
matlab2tikz(fullfile(outDir,'usabilidade_decomposicao.tex'), 'width','\figurewidth','height','\figureheight');

%% ---- Resumo numerico ----
fprintf('\n===== USABILIDADE (perda de SNR em dB; menor=melhor) =====\n');
fprintf('Referencias por SNR:  semComp / manifold(b_hat) / otimo\n');
for iSNR = 1:nSNR
    fprintf('SNR=%+3d dB:  semComp=%+6.2f  manifold=%+6.2f  otimo=%+6.2f\n', ...
        range_SNR_dB(iSNR), loss_nc(iSNR), loss_man(iSNR), loss_opt(iSNR));
    fprintf('            desacopla: ');
    for im=1:nMethods, fprintf('%s=%+6.2f  ', methods{im}, loss_dec(im,iSNR)); end
    fprintf('ORACLE=%+6.2f\n', loss_dec_o(iSNR));
    fprintf('            (ruido dB: ');
    for im=1:nMethods, fprintf('%s=%+5.1f ', methods{im}, 10*log10(mnamp(im,iSNR))); end
    fprintf('| fid: ');
    for im=1:nMethods, fprintf('%s=%.2f ', methods{im}, mfid(im,iSNR)); end
    fprintf(')\n');
end
fprintf(['\nLeitura: desacopla_X < semComp -> compensar AJUDA. ' ...
         'desacopla_X > semComp -> PIORA. Se manifold(b_hat) <= melhor ' ...
         'desacopla -> nao adianta estimar C (use o manifold).\n']);
fprintf('\nFiguras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function [Gdec, namp, fid] = usability(C_hat, C_true, a)
% Ganho de array com desacoplamento D=inv(C_hat) e peso w=a (||a||=1):
%   Gdec = |a^H D C_true a|^2 / (a^H D D^H a)
%   namp = trace(D D^H)/M           (amplificacao de ruido)
%   fid  = |a^H DCa|^2/(||a||^2 ||DCa||^2)   (fidelidade de steering, 0..1)
    M = numel(a);
    if ~all(isfinite(C_hat(:))), Gdec=NaN; namp=NaN; fid=NaN; return; end
    D   = safe_inv(C_hat);
    DCa = D*(C_true*a);
    sig = a'*DCa;
    noi = real(a'*(D*D')*a);
    nDCa= real(DCa'*DCa);
    if noi>eps && isfinite(noi), Gdec = abs(sig)^2/noi; else, Gdec = NaN; end
    if nDCa>eps,                 fid  = abs(sig)^2/(real(a'*a)*nDCa); else, fid = NaN; end
    namp = real(trace(D*D'))/M;
    if ~isfinite(namp), namp = NaN; end
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
        if strcmpi(method,'ORACLE')
            phi = phi_true;
        else
            phi = doa_estimate(Y, method, q, r, lambda, beta_uca, A_dict, phi_grid_deg);
        end
        a_hat = utils.steering_vec_uca(M, r, lambda, theta_deg, phi);
        C_ls  = estimate_C_circulant_uca(b_hat_orig, a_hat, M);
        C_hat = (1-u)*C_hat + u*C_ls;
        frob_it(it) = norm(C_hat - C_true,'fro')/normCtrue;
        cres_it(it) = comp_residual(C_hat, C_true);
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

function rr = comp_residual(C_hat, C_true)
    if ~all(isfinite(C_hat(:))), rr = NaN; return; end
    M = size(C_hat,1);
    E = safe_inv(C_hat) * C_true;
    g = trace(E)/M; if abs(g) > eps, E = E / g; end
    if all(isfinite(E(:))), rr = norm(E - eye(M),'fro')/sqrt(M); else, rr = NaN; end
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
