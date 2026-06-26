% =========================================================================
% DoA_KW_UCA_SelfCal_Iterativo
%
% Self-calibration iterativa por COMPENSACAO de acoplamento mutuo, SEM
% interferente, UMA UNICA direcao por realizacao (single-shot). Estuda a
% CONVERGENCIA do laco alternante:
%
%   Ideia-chave (laco alternante, sem oraculo):
%     1) estima a DoA do sinal a partir do sinal compensado
%     2) estima a matriz de acoplamento C que faz a forma de onda CONHECIDA
%        ter exatamente aquela DoA estimada  (LS circulante)
%     3) compensa (D = inv(C)) e repete
%
%   4 estimadores de DoA dirigindo o mesmo laco, em paralelo:
%       KW (known waveform), Delay-and-Sum, Capon (MPDR) e MUSIC.
%   Todos com init mode = KW.
%
%   MONTE CARLO: varre-se MUITOS azimutes aleatorios (todos no grid de scan)
%   e, em cada azimute, repete-se a medida com varias realizacoes
%   independentes de ruido/dados. Metricas agregadas:
%       - DoA  -> RMSE e MEDIANA (graus) sobre angulos x realizacoes
%       - C    -> erro de Frobenius MEDIO sobre angulos x realizacoes
%
%   Alem das figuras de convergencia (vs iteracao) e do resumo (vs SNR),
%   gera SCATTER de erro POR AZIMUTE: revela se as falhas dos metodos
%   espectrais sao localizadas em certos angulos (e nao uniformes), e
%   compara, por metodo, iterativo vs sem compensacao em cada azimute.
%
% Deliberadamente SEM: criterio de parada, regularizacao de passo, damping,
% adivinhacao/uso da matriz original (C_true so' entra para medir o erro).
% =========================================================================

clear; clc; close all;

%% ---- Parametros do array ----
M      = 8;            % nº de elementos do UCA
fc     = 500e6;        % Hz
c      = 3e8;
lambda = c/fc;
r      = 0.25*lambda;   % raio 0.2 lambda

theta_sig_deg = 90;    % elevacao (plano XY)

%% ---- Parametros do sinal (FSK-2 conhecido) ----
fs    = 288000;        % taxa de amostragem (Hz)
N     = 9900;          % nº de amostras (multiplo de sps)
Rs    = 9600;          % taxa de simbolos
sps   = 30;            % amostras por simbolo
alpha = 0.3;           % roll-off RRC
span  = 8;             % comprimento RRC (simbolos)
fd    = 4.8e3;         % desvio de frequencia FSK

%% ---- Parametros do experimento ----
maxIter      = 2;                  % <-- LIMITE de iteracoes (parametro)
range_SNR_dB = [-6, 0, 6, 12];      % cenarios de SNR
methods      = {'KW','DAS','CAPON','MUSIC'};   % rodam em paralelo
nMethods     = numel(methods);
phi_grid_deg = -180:0.5:180;        % grade do scan p/ DAS/Capon/MUSIC
grid_step    = phi_grid_deg(2) - phi_grid_deg(1);   % passo do grid (0.5 deg)
beta_uca     = 2*pi*(0:M-1).'/M;    % posicoes angulares do UCA (rad)

% --- Monte Carlo ---  (sugestao 3: MUITOS azimutes aleatorios, parametrizado)
n_angles = 36;         % <-- nº de azimutes aleatorios (exercita todas as posicoes)
n_trials = 12;         % <-- nº de medidas (ruido/dados) por azimute
rng(2026, 'twister');  % reprodutibilidade (angulos + ruido)

% Azimutes sorteados em [-180,180), arredondados ao grid (phi fica NO grid).
phi_set = round( (-180 + 360*rand(1, n_angles)) / grid_step ) * grid_step;
fprintf('Azimutes (no grid): %s\n', num2str(sort(phi_set), '%+.1f '));

outDir = fullfile(pwd, 'Iterative SelfCal Graphs');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---- Matriz de acoplamento VERDADEIRA (geometria do UCA) ----
% C_true SO' e' usada para medir o erro de estimacao (nunca dentro do laco).
C_true    = compute_Ctx_for_R(fc, M, r, 50);
normCtrue = norm(C_true, 'fro');
fprintf('C_true gerada (raio=%.2f lambda, ||C_true||_F = %.4f)\n', r/lambda, normCtrue);

%% ---- Dicionario de steering vectors (theta=90) p/ DAS/Capon/MUSIC ----
nGrid  = numel(phi_grid_deg);
A_dict = zeros(M, nGrid);
for ig = 1:nGrid
    A_dict(:,ig) = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_grid_deg(ig));
end

%% ---- Acumuladores ----
nSNR           = numel(range_SNR_dB);
% Convergencia (somas sobre angulos x realizacoes -> RMSE/medio por iteracao)
sumsq_err_iter = zeros(nMethods, maxIter, nSNR);  % SOMA (erro DoA)^2 por iteracao
sum_errF_iter  = zeros(nMethods, maxIter, nSNR);  % SOMA erro Frobenius por iteracao
sum_frob_orc   = zeros(maxIter, nSNR);            % SOMA Frobenius do ORACLE (ang. sabido) por iteracao
cnt            = zeros(1, nSNR);                   % nº de realizacoes por SNR
% Erros por realizacao (final-iter / referencias) p/ scatter, RMSE e MEDIANA
Eabs_final  = zeros(nMethods, n_angles, n_trials, nSNR);  % |erro DoA| iterativo (ultima iter)
Eabs_noComp = zeros(nMethods, n_angles, n_trials, nSNR);  % |erro DoA| sem compensacao
Eabs_noMC   = zeros(nMethods, n_angles, n_trials, nSNR);  % |erro DoA| sem acoplamento

%% ======================= LOOP PRINCIPAL =======================
t_start = tic;
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fprintf('\n===== SNR = %+d dB =====\n', SNR_dB);

    for ia = 1:n_angles
        phi_sig_deg = phi_set(ia);

        for itr = 1:n_trials
            % --- Gera dados; SEM interferente: usamos apenas Xsig e Xn ---
            [~, q, ~, Xsig, ~, Xn] = utils.simulate_fsk_data_uca(M, r, lambda, ...
                phi_sig_deg, phi_sig_deg+90, theta_sig_deg, theta_sig_deg, ...
                SNR_dB, -100, N, fs, Rs, sps, alpha, span, fd);
            q = q(:);

            X_ideal   = Xsig + Xn;          % sem acoplamento  (ruido no receptor)
            X_coupled = C_true*Xsig + Xn;   % com acoplamento  (ruido no receptor)

            % --- DoA de referencia (sem iterar) ---
            for im = 1:nMethods
                phiMC = doa_estimate(X_ideal,   methods{im}, q, r, lambda, beta_uca, A_dict, phi_grid_deg);
                phiNC = doa_estimate(X_coupled, methods{im}, q, r, lambda, beta_uca, A_dict, phi_grid_deg);
                Eabs_noMC(im,ia,itr,iSNR)   = abs(wrapTo180(phiMC - phi_sig_deg));
                Eabs_noComp(im,ia,itr,iSNR) = abs(wrapTo180(phiNC - phi_sig_deg));
            end

            % --- Assinatura espacial acoplada (FIXA ao longo das iteracoes) ---
            qHq        = q'*q;
            b_hat_orig = X_coupled*conj(q)/qHq;

            % --- ORACLE iterativo: MESMO laco, mas usando SEMPRE o angulo
            %     verdadeiro (DoA perfeita). Sem damping (u=1), a cada iteracao
            %     C = LS(b_hat, a(phi_true)) e' constante -> trajetoria plana.
            a_true   = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_sig_deg);
            C_orc    = estimate_C_circulant_uca(b_hat_orig, a_true, M);
            frob_orc = norm(C_orc - C_true,'fro')/normCtrue;
            sum_frob_orc(:,iSNR) = sum_frob_orc(:,iSNR) + frob_orc;   % mesma p/ todas as iteracoes

            % --- Self-cal iterativo: um laco independente por metodo de DoA ---
            for im = 1:nMethods
                mth = methods{im};

                % ===== Init mode = KW =====
                phi0  = doa_estimate(X_coupled, 'KW', q, r, lambda, beta_uca, A_dict, phi_grid_deg);
                a0    = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi0);
                C_hat = estimate_C_circulant_uca(b_hat_orig, a0, M);

                % ===== Laco alternante (DoA -> C -> DoA -> ...) =====
                for it = 1:maxIter
                    D       = inv(C_hat);
                    Y       = D*X_coupled;
                    phi_hat = doa_estimate(Y, mth, q, r, lambda, beta_uca, A_dict, phi_grid_deg);

                    a_hat = utils.steering_vec_uca(M, r, lambda, theta_sig_deg, phi_hat);
                    C_hat = estimate_C_circulant_uca(b_hat_orig, a_hat, M);

                    e_doa = abs(wrapTo180(phi_hat - phi_sig_deg));
                    sumsq_err_iter(im,it,iSNR) = sumsq_err_iter(im,it,iSNR) + e_doa^2;
                    sum_errF_iter(im,it,iSNR)  = sum_errF_iter(im,it,iSNR)  + norm(C_hat - C_true,'fro')/normCtrue;
                    if it == maxIter
                        Eabs_final(im,ia,itr,iSNR) = e_doa;   % erro na ultima iteracao
                    end
                end
            end

            cnt(iSNR) = cnt(iSNR) + 1;
        end
        fprintf('  phi=%+6.1f deg  (%d realizacoes)  | decorrido %.0fs\n', ...
            phi_sig_deg, n_trials, toc(t_start));
    end
end

%% ---- Consolidacao ----
% Convergencia: RMSE de DoA e Frobenius medio por iteracao
rmse_iter  = zeros(nMethods, maxIter, nSNR);
meanF_iter = zeros(nMethods, maxIter, nSNR);
for iSNR = 1:nSNR
    rmse_iter(:,:,iSNR)  = sqrt(sumsq_err_iter(:,:,iSNR) / cnt(iSNR));
    meanF_iter(:,:,iSNR) =      sum_errF_iter(:,:,iSNR)  / cnt(iSNR);
end
mean_frob_orc = sum_frob_orc ./ cnt;   % (maxIter x nSNR): trajetoria do oracle

% Por azimute: RMSE e MEDIANA sobre as realizacoes (dim 3 = trials)
rmse_angle_iter   = squeeze(sqrt(mean(Eabs_final.^2,  3)));  % nMethods x n_angles x nSNR
rmse_angle_noComp = squeeze(sqrt(mean(Eabs_noComp.^2, 3)));
med_angle_iter    = squeeze(median(Eabs_final,  3));
med_angle_noComp  = squeeze(median(Eabs_noComp, 3));

% Global (sobre angulos x realizacoes): RMSE e MEDIANA
E2_final  = reshape(Eabs_final,  nMethods, n_angles*n_trials, nSNR);
E2_noComp = reshape(Eabs_noComp, nMethods, n_angles*n_trials, nSNR);
E2_noMC   = reshape(Eabs_noMC,   nMethods, n_angles*n_trials, nSNR);
rmse_final   = squeeze(rmse_iter(:,maxIter,:));               % nMethods x nSNR
rmse_noComp  = squeeze(sqrt(mean(E2_noComp.^2, 2)));
rmse_noMC    = squeeze(sqrt(mean(E2_noMC.^2,   2)));
med_final    = squeeze(median(E2_final,  2));
med_noComp   = squeeze(median(E2_noComp, 2));
med_noMC     = squeeze(median(E2_noMC,   2));

% Piso para plotagem em eixo log (evita log(0) quando o erro e' nulo)
flr = @(x) max(x, 1e-3);

%% ======================= PLOTS =======================
colorsM   = lines(nMethods);
markers_m = {'o-','s-','d-','^-'};
mk_only   = {'o','s','d','^'};

% --- (A) Uma figura por SNR: RMSE de DoA e Frobenius vs iteracao ---
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fig = figure('Name', sprintf('SelfCal iterativo SNR=%+d dB', SNR_dB), ...
                 'Color','w', 'Position',[80 80 1250 520]);

    subplot(1,2,1); hold on; grid on;
    for im = 1:nMethods
        plot(1:maxIter, flr(squeeze(rmse_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:), 'LineWidth',1.6, 'MarkerFaceColor',colorsM(im,:), ...
            'DisplayName', sprintf('%s (iterativo)', methods{im}));
        yline(flr(rmse_noComp(im,iSNR)), ':', 'Color',colorsM(im,:), ...
            'LineWidth',1.0, 'HandleVisibility','off');
    end
    set(gca,'YScale','log');
    xlabel('Iteracao'); ylabel('RMSE de DoA (graus)');
    title(sprintf('RMSE de DoA vs iteracao  (SNR=%+d dB)', SNR_dB));
    legend('Location','best'); xlim([1 maxIter]);
    subtitle('linha pontilhada (...) = RMSE sem compensacao por metodo');

    subplot(1,2,2); hold on; grid on;
    for im = 1:nMethods
        plot(1:maxIter, flr(squeeze(meanF_iter(im,:,iSNR))), markers_m{im}, ...
            'Color',colorsM(im,:), 'LineWidth',1.6, 'MarkerFaceColor',colorsM(im,:), ...
            'DisplayName',methods{im});
    end
    plot(1:maxIter, flr(mean_frob_orc(:,iSNR)), 'k--p', 'LineWidth',1.6, ...
        'MarkerFaceColor','k','MarkerSize',5, 'DisplayName','oracle (ang. sabido)');
    set(gca,'YScale','log');
    xlabel('Iteracao'); ylabel('media de ||C_{hat}-C_{true}||_F / ||C_{true}||_F');
    title('Erro da matriz de acoplamento vs iteracao');
    legend('Location','best'); xlim([1 maxIter]);

    exportgraphics(fig, fullfile(outDir, ...
        sprintf('selfcal_iter_SNR_%+03d.png', SNR_dB)), 'Resolution',180);
    matlab2tikz(fullfile(outDir, sprintf('selfcal_iter_SNR_%+03d.tex', SNR_dB)), 'width','\figurewidth','height','\figureheight');
end

% --- (B) SCATTER de erro POR AZIMUTE (sugestao 1): 1 figura/SNR, 4 subplots ---
[phi_sorted, ord] = sort(phi_set);
for iSNR = 1:nSNR
    SNR_dB = range_SNR_dB(iSNR);
    fig = figure('Name', sprintf('Erro por azimute SNR=%+d dB', SNR_dB), ...
                 'Color','w', 'Position',[60 60 1300 820]);
    vv   = flr([reshape(rmse_angle_iter(:,:,iSNR),[],1); reshape(rmse_angle_noComp(:,:,iSNR),[],1)]);
    yl_s = [min(vv)*0.8, max(vv)*1.3];   % limites Y COMUNS aos 4 paineis
    for im = 1:nMethods
        subplot(2,2,im); hold on; grid on;
        semilogy(phi_sorted, flr(rmse_angle_iter(im,ord,iSNR)), [mk_only{im} '-'], ...
            'Color',colorsM(im,:), 'LineWidth',1.4, 'MarkerFaceColor',colorsM(im,:), ...
            'DisplayName','iterativo (RMSE/azimute)');
        semilogy(phi_sorted, flr(rmse_angle_noComp(im,ord,iSNR)), [mk_only{im} ':'], ...
            'Color',colorsM(im,:), 'LineWidth',1.0, 'MarkerFaceColor','none', ...
            'DisplayName','sem comp. (RMSE/azimute)');
        yline(flr(med_final(im,iSNR)), 'k--', 'LineWidth',1.0, ...
            'DisplayName','mediana global (iter)');
        set(gca,'YScale','log');
        xlabel('azimute \phi (graus)'); ylabel('RMSE de DoA (graus)');
        xlim([-180 180]); xticks(-180:90:180); ylim(yl_s);
        title(sprintf('%s  | mediana iter=%.3f, semComp=%.3f', ...
            methods{im}, med_final(im,iSNR), med_noComp(im,iSNR)));
        legend('Location','best');
    end
    sgtitle(sprintf('Erro de DoA por azimute  (SNR=%+d dB, %d realizacoes/azimute)', ...
        SNR_dB, n_trials), 'FontWeight','bold');
    exportgraphics(fig, fullfile(outDir, ...
        sprintf('scatter_azimute_SNR_%+03d.png', SNR_dB)), 'Resolution',170);
    matlab2tikz(fullfile(outDir, sprintf('scatter_azimute_SNR_%+03d.tex', SNR_dB)), 'width','\figurewidth','height','\figureheight');
end

% --- (C) Resumo vs SNR: RMSE e MEDIANA (iterativo vs sem compensacao) ---
plot_summary_vs_snr(range_SNR_dB, rmse_final, rmse_noComp, ...
    methods, colorsM, markers_m, flr, n_angles, n_trials, r/lambda, ...
    'RMSE de DoA', fullfile(outDir,'resumo_RMSE_doa_vs_snr.png'));
plot_summary_vs_snr(range_SNR_dB, med_final, med_noComp, ...
    methods, colorsM, markers_m, flr, n_angles, n_trials, r/lambda, ...
    'MEDIANA de erro de DoA', fullfile(outDir,'resumo_MEDIANA_doa_vs_snr.png'));

%% ---- Resumo numerico no terminal ----
fprintf('\n========== RESUMO (erro de DoA, graus) ==========\n');
for iSNR = 1:nSNR
    fprintf('SNR = %+3d dB:\n', range_SNR_dB(iSNR));
    fprintf('   %-6s | %-22s | %-22s | %-12s\n', 'metodo', ...
        'iterativo (RMSE/med)', 'semComp (RMSE/med)', 'semAcop RMSE');
    for im = 1:nMethods
        fprintf('   %-6s | %8.3f / %8.3f   | %8.3f / %8.3f   | %8.3f\n', ...
            methods{im}, rmse_final(im,iSNR), med_final(im,iSNR), ...
            rmse_noComp(im,iSNR), med_noComp(im,iSNR), rmse_noMC(im,iSNR));
    end
end
fprintf('\nFiguras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function plot_summary_vs_snr(range_SNR_dB, M_iter, M_noComp, methods, ...
    colorsM, markers_m, flr, n_angles, n_trials, r_over_lambda, ylab, fname)
% Plota metrica (RMSE ou mediana) vs SNR: iterativo (solido) vs sem comp (pontilhado)
    nMethods = numel(methods);
    fig = figure('Color','w','Position',[100 100 920 560]); hold on; grid on;
    for im = 1:nMethods
        plot(range_SNR_dB, flr(M_iter(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
            'LineWidth',2.0, 'MarkerFaceColor',colorsM(im,:), 'MarkerSize',8, ...
            'DisplayName', sprintf('%s iterativo', methods{im}));
        plot(range_SNR_dB, flr(M_noComp(im,:)), markers_m{im}, 'Color',colorsM(im,:), ...
            'LineWidth',1.2, 'LineStyle',':', 'MarkerFaceColor','none', 'MarkerSize',8, ...
            'DisplayName', sprintf('%s sem comp.', methods{im}));
    end
    set(gca,'YScale','log');
    xlabel('SNR (dB)'); ylabel(sprintf('%s (graus)', ylab));
    xticks(range_SNR_dB); xlim([min(range_SNR_dB)-1, max(range_SNR_dB)+1]);
    title(sprintf(['%s: laco iterativo (solido) vs sem compensacao (pontilhado)\n' ...
                   '%d azimutes aleatorios x %d realizacoes,  raio=%.2f\\lambda'], ...
                   ylab, n_angles, n_trials, r_over_lambda));
    legend('Location','eastoutside');
    exportgraphics(fig, fname, 'Resolution',180);
    matlab2tikz(regexprep(fname,'\.png$','.tex'), 'width','\figurewidth','height','\figureheight');
end

function phi_hat = doa_estimate(X, method, q, r, lambda, beta_uca, A_dict, phi_grid_deg)
% DOA_ESTIMATE  Estima o azimute (graus) por um dos 4 metodos.
    M = size(X,1);
    switch upper(method)
        case 'KW'
            [~, phi] = doa_kw_uca(X, q(:).', r, lambda, beta_uca);
            phi_hat  = phi(1);
        case 'DAS'
            R      = (X*X')/size(X,2);
            scores = real(sum(conj(A_dict).*(R*A_dict), 1));
            [~,ip] = max(scores);
            phi_hat = phi_grid_deg(ip);
        case 'CAPON'
            R      = (X*X')/size(X,2);
            R      = R + 1e-6*trace(R)/M*eye(M);
            Rinv   = R\eye(M);
            den    = real(sum(conj(A_dict).*(Rinv*A_dict), 1));
            scores = 1./max(den, eps);
            [~,ip] = max(scores);
            phi_hat = phi_grid_deg(ip);
        case 'MUSIC'
            R = (X*X')/size(X,2);
            R = (R+R')/2;
            [V, Dg]   = eig(R);
            [~, idx]  = sort(real(diag(Dg)), 'descend');
            V         = V(:, idx);
            P_sources = 1;                 % sem interferente: 1 fonte
            En        = V(:, P_sources+1:end);
            proj      = En'*A_dict;
            den       = real(sum(conj(proj).*proj, 1));
            scores    = 1./max(den, eps);
            [~,ip]    = max(scores);
            phi_hat   = phi_grid_deg(ip);
        otherwise
            error('doa_estimate:method', 'Metodo desconhecido: %s', method);
    end
end

function Ctx = compute_Ctx_for_R(fc, M, R, Z0)
% COMPUTE_CTX_FOR_R  Matriz de acoplamento de um UCA de dipolos (parametros-S).
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
