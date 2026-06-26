% =========================================================================
% DoA_KW_UCA_ISR_Sweep
%
% Demonstra o problema da ESTIMACAO DE DoA com interferente forte:
% metodos classicos (DAS/Capon/MUSIC) sao baseados em POTENCIA/SUBESPACO e
% travam na fonte MAIS FORTE. Quando ISR > 0 (interferente mais forte que o
% SOI), eles apontam para o INTERFERENTE, nao para o SOI.
%
% O KW (known-waveform) usa a forma de onda conhecida do preambulo: a
% correlacao cruzada b_hat = X q*/(q'q) ~ C a_SOI e' INDEPENDENTE da potencia
% do interferente (descorrelacionado de q). Por isso o KW continua travado no
% SOI mesmo com ISR >> 0.
%
% Cenario: 1 SOI (angulo fixo por realizacao) + 1 interferente co-canal a
% sep_deg de distancia, acoplamento mutuo no canal, SNR fixo. Varre-se o ISR.
%
% Metricas vs ISR (por metodo):
%   (1) RMSE do azimute do SOI (graus)   -> dispara quando o metodo trava no interf.
%   (2) P(travar no interferente) (%)    -> fracao de estimativas mais perto do interf.
% =========================================================================

clear; clc; close all;
warning('off','MATLAB:singularMatrix'); warning('off','MATLAB:nearlySingularMatrix');

%% ---- Parametros ----
M = 8; fc = 500e6; c = 3e8; lambda = c/fc; r = 0.2*lambda; theta_sig_deg = 90;
fs = 288000; Rs = 9600; sps = 30; alpha = 0.3; span = 8; fd = 4.8e3;

K_kw     = 2100;            % preambulo (amostras p/ estimar DoA); sem payload aqui
K        = K_kw; k0 = 1;
SNR_dB   = 10;             % SNR do SOI FIXO (o problema e' o ISR, nao o ruido)
sep_deg  = 40;            % separacao SOI-interferente (grande p/ o "trava" ficar visivel)
interf_t = 'fsk2';        % interferente co-canal (mesma modulacao)

range_ISR_dB = [-15 -10 -6 -3 0 3 6 9 12 15 20];
methods = {'KW','DAS','CAPON','MUSIC'};  nM = numel(methods);
phi_grid_deg = -180:0.5:180;
beta_uca = 2*pi*(0:M-1).'/M;

% Monte Carlo
n_angles = 40; n_trials = 10; n_real = n_angles*n_trials;
rng(2026,'twister');
gstep = phi_grid_deg(2)-phi_grid_deg(1);
phi_set = round((-180 + 360*rand(1,n_angles))/gstep)*gstep;

outDir = fullfile(pwd,'ISR Sweep Graphs'); if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---- C_true (acoplamento) ----
C_true = compute_Ctx_for_R(fc, M, r, 50);
fprintf('ISR sweep | raio=%.2f lambda, cond(C)=%.1f, SNR=%+d dB, sep=%d deg, K_kw=%d\n', ...
    r/lambda, cond(C_true), SNR_dB, sep_deg, K_kw);

%% ---- Acumuladores ----
nISR = numel(range_ISR_dB);
sumsq_soi = zeros(nM,nISR);    % soma dos erros^2 ao SOI (p/ RMSE)
cnt_lock  = zeros(nM,nISR);    % nº de vezes travado no interferente
cnt       = zeros(1,nISR);

t0 = tic;
for iISR = 1:nISR
    ISR_dB = range_ISR_dB(iISR);
    for ia = 1:n_angles
        phi_sig = phi_set(ia);
        phi_int = wrapTo180(phi_sig + sep_deg);
        for itr = 1:n_trials
            [Xtot, y_ref, k0t, Xsig, Xint, Xn] = utils.simulate_data_uca_v3(M, r, lambda, ...
                phi_sig, phi_int, theta_sig_deg, theta_sig_deg, SNR_dB, ISR_dB, K, K_kw, ...
                fs, Rs, sps, alpha, span, fd, interf_t, k0);
            q = y_ref(:);  win = k0t : k0t + K_kw - 1;
            Xk_c = (C_true*(Xsig + Xint) + Xn);  Xk_c = Xk_c(:, win);   % preambulo acoplado

            for m = 1:nM
                phi_hat = doa_estimate(Xk_c, methods{m}, q, r, lambda, beta_uca, phi_grid_deg);
                err_soi = abs(wrapTo180(phi_hat - phi_sig));
                err_int = abs(wrapTo180(phi_hat - phi_int));
                sumsq_soi(m,iISR) = sumsq_soi(m,iISR) + err_soi^2;
                if err_int < err_soi, cnt_lock(m,iISR) = cnt_lock(m,iISR) + 1; end
            end
            cnt(iISR) = cnt(iISR) + 1;
        end
    end
    fprintf('  ISR=%+3d dB feito (%.0fs)\n', ISR_dB, toc(t0));
end

%% ---- Consolidacao ----
RMSE_soi = sqrt(sumsq_soi ./ cnt);     % nM x nISR
Plock    = 100 * cnt_lock ./ cnt;      % %

colorsM = lines(nM);  markers_m = {'o-','s-','d-','^-'};

%% =================== FIGURA: RMSE e P(trava) vs ISR ===================
fig = figure('Color','w','Position',[50 60 1300 520]);

subplot(1,2,1); hold on; grid on;
for m = 1:nM
    plot(range_ISR_dB, RMSE_soi(m,:), markers_m{m}, 'Color',colorsM(m,:), ...
        'LineWidth',1.9,'MarkerFaceColor',colorsM(m,:),'MarkerSize',7,'DisplayName',methods{m});
end
yline(sep_deg,'k:','LineWidth',1.4,'DisplayName','erro \approx sep (travou no interf.)');
xline(0,'Color',[.5 .5 .5],'LineWidth',1.0,'HandleVisibility','off');   % ISR=0
set(gca,'YScale','log'); xlabel('ISR (dB)  [interferente - SOI]'); ylabel('RMSE do azimute do SOI (graus)');
xticks(range_ISR_dB); title('Erro de apontamento ao SOI vs ISR'); legend('Location','northwest');

subplot(1,2,2); hold on; grid on;
for m = 1:nM
    plot(range_ISR_dB, Plock(m,:), markers_m{m}, 'Color',colorsM(m,:), ...
        'LineWidth',1.9,'MarkerFaceColor',colorsM(m,:),'MarkerSize',7,'DisplayName',methods{m});
end
yline(50,'k:','LineWidth',1.0,'HandleVisibility','off');
xline(0,'Color',[.5 .5 .5],'LineWidth',1.0,'HandleVisibility','off');
xlabel('ISR (dB)  [interferente - SOI]'); ylabel('P(travar no interferente) (%)');
xticks(range_ISR_dB); ylim([-3 103]); title('Fracao de estimativas que travam no interferente'); legend('Location','northwest');

sgtitle(sprintf(['DoA com interferente forte: classicos travam na fonte mais forte; KW (known-waveform) nao\n' ...
    'raio=%.2f\\lambda, SNR=%+d dB, sep=%d^o, %d azimutes x %d realizacoes'], ...
    r/lambda, SNR_dB, sep_deg, n_angles, n_trials),'FontWeight','bold');
exportgraphics(fig, fullfile(outDir,'isr_sweep_doa.png'),'Resolution',170);
matlab2tikz(fullfile(outDir,'isr_sweep_doa.tex'), 'width','\figurewidth','height','\figureheight');

%% ---- Resumo numerico ----
fprintf('\n===== P(travar no interferente) %% vs ISR =====\n');
fprintf('%-8s', 'ISR'); fprintf('%8d', range_ISR_dB); fprintf('\n');
for m = 1:nM, fprintf('%-8s', methods{m}); fprintf('%8.1f', Plock(m,:)); fprintf('\n'); end
fprintf('\n===== RMSE do SOI (graus) vs ISR =====\n');
fprintf('%-8s', 'ISR'); fprintf('%8d', range_ISR_dB); fprintf('\n');
for m = 1:nM, fprintf('%-8s', methods{m}); fprintf('%8.2f', RMSE_soi(m,:)); fprintf('\n'); end
fprintf('\nFiguras salvas em: %s\n', outDir);

% =========================================================================
%                            FUNCOES LOCAIS
% =========================================================================
function phi_hat = doa_estimate(X, method, q, r, lambda, beta_uca, phi_grid_deg)
% DoA do SOI. KW usa a forma de onda conhecida (robusto ao interferente);
% DAS/Capon/MUSIC varrem o espectro assumindo 1 FONTE dominante (a mais forte).
    M = size(X,1);
    persistent A_dict_cache phi_cache r_cache
    if isempty(phi_cache) || ~isequal(phi_cache,phi_grid_deg) || ~isequal(r_cache,r)
        A_dict_cache = zeros(M,numel(phi_grid_deg));
        for ig=1:numel(phi_grid_deg), A_dict_cache(:,ig)=utils.steering_vec_uca(M,r,lambda,90,phi_grid_deg(ig)); end
        phi_cache = phi_grid_deg;  r_cache = r;
    end
    A_dict = A_dict_cache;
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
            En = V(:, 2:end);            % 1 FONTE (estima so' a fonte dominante)
            proj = En'*A_dict;
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
