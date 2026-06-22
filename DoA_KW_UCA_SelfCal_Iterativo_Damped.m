% =========================================================================
% DoA_KW_UCA_SelfCal_Iterativo_Damped
%
% Derivado de DoA_KW_UCA_SelfCal_Iterativo (single-shot, uma direcao por
% realizacao). Acrescenta:
%
%   (1) PASSO DE RELAXACAO 'u' na evolucao de C (suavizacao / sub-relaxacao):
%         C_t = (1-u)*C_{t-1} + u*C_LS
%       u = 1  -> comportamento original (sem suavizacao)
%       u < 1  -> suaviza a trajetoria de C, atenua oscilacoes/ciclos-limite
%
%   (2) METRICAS de C alem do Frobenius bruto (que e' enganoso):
%         - Frobenius invariante a escala
%         - Residuo de compensacao efetiva  ||norm(inv(C_hat) C_true) - I||/sqrt(M)
%         - Erro por coeficiente |c_hat_k - c_k| e PLANO COMPLEXO (verdadeiro
%           vs nuvem de estimativas), K = floor(M/2) coeficientes unicos
%
%   (3) REFERENCIA "oracle de angulo conhecido": estima C por LS circulante
%       usando a DIRECAO VERDADEIRA a(phi_true) (1 shot, sem iterar). E' o
%       melhor C alcancavel dado o ruido de b_hat e o piso estrutural da
%       geometria -> limite inferior do erro de C.
%
% Mantem 4 metodos de DoA em paralelo (KW/DAS/Capon/MUSIC), init = KW.
% =========================================================================

clear; clc; close all;

%% ---- Parametros do array ----
M      = 8;            % nº de elementos do UCA
fc     = 500e6;
c      = 3e8;
lambda = c/fc;
r      = 0.1*lambda;   % raio (abertura pequena -> acoplamento forte, C dificil)
theta_sig_deg = 90;

%% ---- Parametros do sinal (FSK-2 conhecido) ----
fs = 288000; N = 9900;%2100; 
Rs = 9600; sps = 30; alpha = 0.3; span = 8; fd = 4.8e3;

%% ---- Parametros do experimento ----
maxIter      = 2;
u            = 1;                 % <-- PASSO de relaxacao da atualizacao de C (1 = sem suavizar)
range_u      = [1.0 0.7 0.5 0.3];   % varredura de u (figura dedicada)
snr_u_demo   = 6;                   % SNR usada na varredura de u (dB)
n_real_u     = 150;                 % realizacoes na varredura de u

range_SNR_dB = [-6, 0, 6, 12];
methods      = {'KW','DAS','CAPON','MUSIC'};
nMethods     = numel(methods);
phi_grid_deg = -180:0.5:180;
grid_step    = phi_grid_deg(2) - phi_grid_deg(1);
beta_uca     = 2*pi*(0:M-1).'/M;
K            = floor(M/2);

% --- Monte Carlo ---
n_angles = 24;        % azimutes aleatorios (no grid)
n_trials = 10;        % realizacoes por azimute
n_real   = n_angles*n_trials;
rng(2026, 'twister');
phi_set  = round( (-180 + 360*rand(1, n_angles)) / grid_step ) * grid_step;

outDir = fullfile(pwd, 'Iterative Damped Graphs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---- Matriz de acoplamento VERDADEIRA (so' p/ medir erro) ----
C_true     = compute_Ctx_for_R(fc, M, r, 50);
normCtrue  = norm(C_true, 'fro');
c_true_vec = C_true(1, 2:K+1).';
fprintf('C_true (raio=%.2f lambda, ||C_true||_F=%.4f, cond=%.1f),  u=%.2f\n', ...
    r/lambda, normCtrue, cond(C_true), u);

%% ---- Dicionario de steering vectors ----
nGrid = numel(phi_grid_deg);
A_dict = zeros(M, nGrid);
for ig = 1:nGrid
    A_dict(:,ig) = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_grid_deg(ig));
end

%% =====================================================================
%   PARTE A: convergencia, metricas de C e oracle (passo u fixo)
%% =====================================================================
nSNR             = numel(range_SNR_dB);
sumsq_doa_iter   = zeros(nMethods, maxIter, nSNR);
sum_frob_iter    = zeros(nMethods, maxIter, nSNR);
sum_cres_iter    = zeros(nMethods, maxIter, nSNR);
sumsq_doa_noComp = zeros(nMethods, nSNR);
c_runs           = zeros(nMethods, K, n_real, nSNR);
% Oracle (angulo sabido) -> agora TRAJETORIA por iteracao (mesma do metodo)
sum_frob_oracle  = zeros(maxIter, nSNR);
sum_cres_oracle  = zeros(maxIter, nSNR);
c_oracle         = zeros(K, n_real, nSNR);
cnt              = zeros(1, nSNR);

t0 = tic;
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fprintf('\n===== SNR = %+d dB  (u=%.2f) =====\n', SNR_dB, u);
    ir = 0;
    for ia = 1:n_angles
        phi_true = phi_set(ia);
        for itr = 1:n_trials
            ir = ir + 1;
            [~, q, ~, Xsig, ~, Xn] = utils.simulate_fsk_data_uca(M, r, lambda, ...
                phi_true, phi_true+90, theta_sig_deg, theta_sig_deg, ...
                SNR_dB, -100, N, fs, Rs, sps, alpha, span, fd);
            q = q(:);
            X_coupled = C_true*Xsig + Xn;
            b_hat_orig = X_coupled*conj(q)/(q'*q);

            % --- DoA de referencia sem compensacao ---
            for im = 1:nMethods
                phiNC = doa_estimate(X_coupled, methods{im}, q, r, lambda, beta_uca, A_dict, phi_grid_deg);
                sumsq_doa_noComp(im,iSNR) = sumsq_doa_noComp(im,iSNR) + wrapTo180(phiNC - phi_true)^2;
            end

            % --- ORACLE iterativo: MESMO laco dos metodos (init KW, passo u),
            %     mas usando SEMPRE o angulo verdadeiro (DoA perfeita) no lugar
            %     da DoA estimada. Vira uma trajetoria comparavel as dos metodos.
            [frob_orc_it, cres_orc_it, ~, c_orc_fin] = selfcal_one_damped('ORACLE', ...
                X_coupled, b_hat_orig, q, phi_true, M, r, lambda, theta_sig_deg, ...
                beta_uca, A_dict, phi_grid_deg, C_true, normCtrue, u, maxIter);
            sum_frob_oracle(:,iSNR) = sum_frob_oracle(:,iSNR) + frob_orc_it(:);
            sum_cres_oracle(:,iSNR) = sum_cres_oracle(:,iSNR) + cres_orc_it(:);
            c_oracle(:,ir,iSNR)     = c_orc_fin(:);

            % --- Self-cal iterativo DAMPED, por metodo ---
            for im = 1:nMethods
                [frob_it, cres_it, derr_it, c_fin] = selfcal_one_damped(methods{im}, ...
                    X_coupled, b_hat_orig, q, phi_true, M, r, lambda, theta_sig_deg, ...
                    beta_uca, A_dict, phi_grid_deg, C_true, normCtrue, u, maxIter);
                sum_frob_iter(im,:,iSNR)  = sum_frob_iter(im,:,iSNR)  + frob_it(:).';
                sum_cres_iter(im,:,iSNR)  = sum_cres_iter(im,:,iSNR)  + cres_it(:).';
                sumsq_doa_iter(im,:,iSNR) = sumsq_doa_iter(im,:,iSNR) + (derr_it(:).').^2;
                c_runs(im,:,ir,iSNR)      = c_fin(:).';
            end
            cnt(iSNR) = cnt(iSNR) + 1;
        end
        fprintf('  phi=%+6.1f deg  | decorrido %.0fs\n', phi_true, toc(t0));
    end
end

% --- Consolidacao ---
rmse_doa_iter  = zeros(nMethods, maxIter, nSNR);
mean_frob_iter = zeros(nMethods, maxIter, nSNR);
mean_cres_iter = zeros(nMethods, maxIter, nSNR);
for iSNR = 1:nSNR
    rmse_doa_iter(:,:,iSNR)  = sqrt(sumsq_doa_iter(:,:,iSNR) / cnt(iSNR));
    mean_frob_iter(:,:,iSNR) =      sum_frob_iter(:,:,iSNR)  / cnt(iSNR);
    mean_cres_iter(:,:,iSNR) =      sum_cres_iter(:,:,iSNR)  / cnt(iSNR);
end
rmse_doa_final = squeeze(rmse_doa_iter(:,maxIter,:));
frob_final     = squeeze(mean_frob_iter(:,maxIter,:));
cres_final     = squeeze(mean_cres_iter(:,maxIter,:));
rmse_doa_noComp= sqrt(sumsq_doa_noComp ./ cnt);
mean_frob_oracle = sum_frob_oracle ./ cnt;   % (maxIter x nSNR): trajetoria do oracle
mean_cres_oracle = sum_cres_oracle ./ cnt;
frob_oracle    = mean_frob_oracle(end,:);    % valor final (p/ resumo vs SNR e tabela)
cres_oracle    = mean_cres_oracle(end,:);

% Metricas derivadas dos coeficientes (Frobenius inv. a escala + erro por coef.)
scaleinv_final = zeros(nMethods, nSNR);
coeff_err      = zeros(nMethods, K, nSNR);
scaleinv_orc   = zeros(1, nSNR);
coeff_err_orc  = zeros(K, nSNR);
for iSNR = 1:nSNR
    accO = 0;
    for run = 1:n_real
        Corc = reconstruct_C(c_oracle(:,run,iSNR), M);
        accO = accO + scaleinv_frob(Corc, C_true);
    end
    scaleinv_orc(iSNR) = accO / n_real;
    for k = 1:K
        coeff_err_orc(k,iSNR) = mean(abs(c_oracle(k,:,iSNR).' - c_true_vec(k)));
    end
    for im = 1:nMethods
        acc = 0;
        for run = 1:n_real
            Chat = reconstruct_C(squeeze(c_runs(im,:,run,iSNR)).', M);
            acc  = acc + scaleinv_frob(Chat, C_true);
        end
        scaleinv_final(im,iSNR) = acc / n_real;
        for k = 1:K
            coeff_err(im,k,iSNR) = mean(abs(squeeze(c_runs(im,k,:,iSNR)) - c_true_vec(k)));
        end
    end
end

flr = @(x) max(x, 1e-3);
colorsM   = lines(nMethods);
markers_m = {'o-','s-','d-','^-'};

%% =====================================================================
%   PARTE U: varredura do passo de relaxacao u  (SNR = snr_u_demo)
%% =====================================================================
fprintf('\n===== [U] Varredura de u em SNR=%+d dB =====\n', snr_u_demo);
nU = numel(range_u);
sum_frob_u = zeros(nMethods, maxIter, nU);
cnt_u = 0;
for run = 1:n_real_u
    phi_true = round((-180 + 360*rand)/grid_step)*grid_step;
    [~, q, ~, Xsig, ~, Xn] = utils.simulate_fsk_data_uca(M, r, lambda, ...
        phi_true, phi_true+90, theta_sig_deg, theta_sig_deg, ...
        snr_u_demo, -100, N, fs, Rs, sps, alpha, span, fd);
    q = q(:); X_coupled = C_true*Xsig + Xn; b_hat_orig = X_coupled*conj(q)/(q'*q);
    for iu = 1:nU
        for im = 1:nMethods
            [frob_it, ~, ~, ~] = selfcal_one_damped(methods{im}, X_coupled, b_hat_orig, ...
                q, phi_true, M, r, lambda, theta_sig_deg, beta_uca, A_dict, phi_grid_deg, ...
                C_true, normCtrue, range_u(iu), maxIter);
            sum_frob_u(im,:,iu) = sum_frob_u(im,:,iu) + frob_it(:).';
        end
    end
    cnt_u = cnt_u + 1;
end
mean_frob_u = sum_frob_u / cnt_u;

%% ======================= PLOTS =======================
% --- (A) convergencia por SNR: DoA, Frobenius+compRes, com oracle ---
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fig = figure('Name',sprintf('Damped SNR=%+d dB',SNR_dB),'Color','w','Position',[40 60 1500 470]);

    subplot(1,3,1); hold on; grid on;
    for im = 1:nMethods
        plot(1:maxIter, flr(squeeze(rmse_doa_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:),'LineWidth',1.6,'MarkerFaceColor',colorsM(im,:),'DisplayName',methods{im});
        yline(flr(rmse_doa_noComp(im,iSNR)), ':','Color',colorsM(im,:),'LineWidth',1.6,'HandleVisibility','off');
    end
    set(gca,'YScale','log'); xlabel('Iteracao'); ylabel('RMSE de DoA (graus)');
    title('DoA (pontilhado = sem comp.)'); legend('Location','best'); xlim([1 maxIter]);

    subplot(1,3,2); hold on; grid on;
    for im = 1:nMethods
        plot(1:maxIter, flr(squeeze(mean_frob_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:),'LineWidth',1.6,'MarkerFaceColor',colorsM(im,:),'DisplayName',methods{im});
    end
    plot(1:maxIter, flr(mean_frob_oracle(:,iSNR)), 'k--p', 'LineWidth',1.6, ...
        'MarkerFaceColor','k','MarkerSize',5, 'DisplayName','oracle (ang. sabido)');
    set(gca,'YScale','log'); xlabel('Iteracao'); ylabel('Frobenius ||C_{hat}-C_{true}||/||C_{true}||');
    title('Frobenius de C vs iteracao'); legend('Location','best'); xlim([1 maxIter]);

    subplot(1,3,3); hold on; grid on;
    for im = 1:nMethods
        plot(1:maxIter, flr(squeeze(mean_cres_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:),'LineWidth',1.6,'MarkerFaceColor',colorsM(im,:),'DisplayName',methods{im});
    end
    plot(1:maxIter, flr(mean_cres_oracle(:,iSNR)), 'k--p', 'LineWidth',1.6, ...
        'MarkerFaceColor','k','MarkerSize',5, 'DisplayName','oracle (ang. sabido)');
    set(gca,'YScale','log'); xlabel('Iteracao'); ylabel('residuo de compensacao');
    title('Residuo de compensacao vs iteracao'); legend('Location','best'); xlim([1 maxIter]);

    sgtitle(sprintf('Damped self-cal  SNR=%+d dB, u=%.2f, raio=%.2f\\lambda', SNR_dB, u, r/lambda),'FontWeight','bold');
    exportgraphics(fig, fullfile(outDir, sprintf('damped_iter_SNR_%+03d.png',SNR_dB)),'Resolution',170);
end

% --- (U) efeito do passo u na trajetoria de C (Frobenius vs iter) ---
u_colors = cool(nU);
vu = flr(mean_frob_u(:));                          % limites Y COMUNS aos 4 paineis
yl_u = [min(vu)*0.9, max(vu)*1.15];
fig = figure('Name','Efeito de u','Color','w','Position',[60 60 1150 820]);
for im = 1:nMethods
    subplot(2,2,im); hold on; grid on;
    for iu = 1:nU
        plot(1:maxIter, flr(squeeze(mean_frob_u(im,:,iu))), 'o-', 'Color',u_colors(iu,:), ...
            'LineWidth',1.6,'MarkerFaceColor',u_colors(iu,:),'DisplayName',sprintf('u=%.1f',range_u(iu)));
    end
    set(gca,'YScale','log'); xlabel('Iteracao'); ylabel('Frobenius de C');
    title(methods{im}); legend('Location','best'); xlim([1 maxIter]); ylim(yl_u);
end
sgtitle(sprintf('Efeito do passo de relaxacao u na evolucao de C  (SNR=%+d dB, raio=%.2f\\lambda)', ...
    snr_u_demo, r/lambda),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'efeito_passo_u.png'),'Resolution',170);

% --- (D) resumo metricas de C vs SNR (3 metricas), com oracle ---
fig = figure('Color','w','Position',[60 80 1500 470]);
metsD = {frob_final, scaleinv_final, cres_final};
orcD  = {frob_oracle, scaleinv_orc, cres_oracle};
titD  = {'Frobenius bruto', 'Frobenius invariante a escala', 'Residuo de compensacao'};
for sp = 1:3
    subplot(1,3,sp); hold on; grid on;
    for im = 1:nMethods
        plot(range_SNR_dB, flr(metsD{sp}(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
            'LineWidth',1.8,'MarkerFaceColor',colorsM(im,:),'MarkerSize',8,'DisplayName',methods{im});
    end
    plot(range_SNR_dB, flr(orcD{sp}), 'k--p', 'LineWidth',1.6,'MarkerFaceColor','k','MarkerSize',8, ...
        'DisplayName','oracle (ang. sabido)');
    set(gca,'YScale','log'); xlabel('SNR (dB)'); ylabel('erro de C'); xticks(range_SNR_dB);
    title(titD{sp}); legend('Location','best');
end
sgtitle(sprintf('Metricas de C vs SNR  (u=%.2f, %d realizacoes, raio=%.2f\\lambda)', ...
    u, n_real, r/lambda),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'resumo_metricas_C_vs_snr.png'),'Resolution',170);

% --- (E) plano complexo dos coeficientes por SNR (metodos + oracle) ---
ck_colors = [0.85 0.10 0.10; 0.10 0.60 0.15; 0.10 0.30 0.85; 0.60 0.15 0.65];
for iSNR = 1:nSNR
    fig = figure('Name',sprintf('Coef. complexos SNR=%+d',range_SNR_dB(iSNR)),'Color','w','Position',[40 40 1500 820]);
    srcC = cell(1, nMethods+1);
    for im = 1:nMethods, srcC{im} = squeeze(c_runs(im,:,:,iSNR)); end
    srcC{nMethods+1} = squeeze(c_oracle(:,:,iSNR));
    labC = [methods, {'ORACLE (ang. sabido)'}];
    errC = [coeff_err(:,:,iSNR); coeff_err_orc(:,iSNR).'];   % (nMethods+1) x K
    lim  = common_complex_lim(srcC, c_true_vec);             % limites COMUNS aos paineis
    for pp = 1:(nMethods+1)
        subplot(2,3,pp); hold on; grid on; axis equal;
        S = srcC{pp};
        for k = 1:K
            scatter(real(S(k,:)), imag(S(k,:)), 16, ck_colors(k,:), 'filled', ...
                'MarkerFaceAlpha',0.30, 'DisplayName',sprintf('c_%d (est)',k));
        end
        for k = 1:K
            plot(real(c_true_vec(k)), imag(c_true_vec(k)), 'p', 'MarkerSize',16, ...
                'MarkerFaceColor',ck_colors(k,:),'MarkerEdgeColor','k','LineWidth',1.3,'HandleVisibility','off');
            text(real(c_true_vec(k)), imag(c_true_vec(k)), sprintf(' c_%d',k),'FontWeight','bold','FontSize',10);
        end
        xlabel('Re'); ylabel('Im'); xlim(lim); ylim(lim);
        estr = strjoin(arrayfun(@(k) sprintf('%.2f',errC(pp,k)), 1:K, 'uni',0), ', ');
        title(sprintf('%s  |\\Deltac_k|=[%s]', labC{pp}, estr), 'FontSize',9);
        if pp == 1, legend('Location','bestoutside'); end
    end
    sgtitle(sprintf(['Coeficientes c_1..c_%d: verdadeiro (estrela) vs estimado (nuvem)\n' ...
        'SNR=%+d dB, u=%.2f, %d realizacoes, raio=%.2f\\lambda'], ...
        K, range_SNR_dB(iSNR), u, n_real, r/lambda),'FontWeight','bold');
    exportgraphics(fig, fullfile(outDir, sprintf('coeficientes_complexo_SNR_%+03d.png',range_SNR_dB(iSNR))),'Resolution',160);
end

%% ---- Resumo numerico ----
fprintf('\n===== RESUMO (u=%.2f, %d realizacoes) =====\n', u, n_real);
for iSNR = 1:nSNR
    fprintf('SNR=%+3d dB:   [oracle] Frob=%.3f FrobEsc=%.3f compRes=%.3f\n', ...
        range_SNR_dB(iSNR), frob_oracle(iSNR), scaleinv_orc(iSNR), cres_oracle(iSNR));
    for im = 1:nMethods
        fprintf('   %-6s  Frob=%.3f FrobEsc=%.3f compRes=%.3f  DoA(rmse)=%.3f (semComp=%.3f)\n', ...
            methods{im}, frob_final(im,iSNR), scaleinv_final(im,iSNR), cres_final(im,iSNR), ...
            rmse_doa_final(im,iSNR), rmse_doa_noComp(im,iSNR));
    end
end
fprintf('\nFiguras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function [frob_it, cres_it, derr_it, c_fin] = selfcal_one_damped(method, X_coupled, ...
    b_hat_orig, q, phi_true, M, r, lambda, theta_deg, beta_uca, A_dict, phi_grid_deg, ...
    C_true, normCtrue, u, maxIter)
% SELFCAL_ONE_DAMPED  Laco alternante single-shot com passo de relaxacao u:
%   C_t = (1-u)*C_{t-1} + u*C_LS.  init = KW.
    % init KW
    phi0  = doa_estimate(X_coupled, 'KW', q, r, lambda, beta_uca, A_dict, phi_grid_deg);
    a0    = utils.steering_vec_uca(M, r, lambda, theta_deg, phi0);
    C_hat = estimate_C_circulant_uca(b_hat_orig, a0, M);

    frob_it=zeros(maxIter,1); cres_it=zeros(maxIter,1); derr_it=zeros(maxIter,1);
    for it = 1:maxIter
        D   = safe_inv(C_hat);
        Y   = D * X_coupled;
        if strcmpi(method,'ORACLE')
            phi = phi_true;     % DoA PERFEITA: sempre o angulo verdadeiro
        else
            phi = doa_estimate(Y, method, q, r, lambda, beta_uca, A_dict, phi_grid_deg);
        end
        a_hat = utils.steering_vec_uca(M, r, lambda, theta_deg, phi);
        C_ls  = estimate_C_circulant_uca(b_hat_orig, a_hat, M);

        % --- passo de relaxacao (suavizacao) ---
        C_hat = (1-u)*C_hat + u*C_ls;

        frob_it(it) = norm(C_hat - C_true,'fro')/normCtrue;
        cres_it(it) = comp_residual(C_hat, C_true);
        derr_it(it) = abs(wrapTo180(phi - phi_true));
    end
    c_fin = C_hat(1, 2:floor(M/2)+1).';
end

function lim = common_complex_lim(srcC, c_true_vec)
% Limites Re/Im COMUNS (quadrados) p/ os paineis do plano complexo, robustos
% a outliers (percentis 2-98%, garantindo que as estrelas verdadeiras caibam).
    parts = cell(numel(srcC)+1,1);
    for pp = 1:numel(srcC)
        S = srcC{pp};
        parts{pp} = [real(S(:)); imag(S(:))];
    end
    parts{end} = [real(c_true_vec(:)); imag(c_true_vec(:))];
    v = vertcat(parts{:});  v = sort(v(isfinite(v)));  n = numel(v);
    if n < 5, lim = [-1 1]; return; end
    lo = v(max(1,round(0.02*n)));  hi = v(min(n,round(0.98*n)));
    lo = min(lo, min(parts{end}));  hi = max(hi, max(parts{end}));
    pad = 0.1*(hi-lo) + eps;  lim = [lo-pad, hi+pad];
end

function C = reconstruct_C(c, M)
    K = floor(M/2); is_even = (mod(M,2)==0); c = c(:);
    if is_even, first_row = [1, c(1:K-1).', c(K), flip(c(1:K-1).')];
    else,       first_row = [1, c(1:K).', flip(c(1:K).')]; end
    C = zeros(M);
    for i = 1:M, C(i,:) = circshift(first_row, [0, i-1]); end
end

function rr = comp_residual(C_hat, C_true)
% Residuo de compensacao efetiva (escala removida). 0=perfeito, ~1=sem compensar.
    M = size(C_hat,1);
    E = safe_inv(C_hat) * C_true;
    g = trace(E)/M; if abs(g) > eps, E = E / g; end
    if all(isfinite(E(:))), rr = min(norm(E - eye(M),'fro')/sqrt(M), 100);
    else,                   rr = 100; end
end

function e = scaleinv_frob(C_hat, C_true)
    s = (C_true(:)' * C_hat(:)) / (C_true(:)' * C_true(:));
    e = norm(C_hat - s*C_true,'fro') / norm(C_true,'fro');
end

function D = safe_inv(C)
    M = size(C,1);
    ws = warning('off','MATLAB:singularMatrix');
    wn = warning('off','MATLAB:nearlySingularMatrix');
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
