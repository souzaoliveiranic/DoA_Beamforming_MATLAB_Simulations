% =========================================================================
% BLOCO: Auto-calibracao CEGA multi-transmissao (KW) — modelo real do script
% -------------------------------------------------------------------------
% COLE este bloco DENTRO do script principal, APOS o bloco de calibracao
% (depois da linha que define C_true = Coupling_matrices_real(:,:,1) e onde
%  q_cal ja' existe). Ele reusa as variaveis ja' carregadas:
%   M, lambda, fs, N, Rs, sps, alpha, span, fd, theta_cal_deg, range_radius,
%   Coupling_matrices_real, outDir, e o preambulo de calibracao q_cal.
%
% Estuda recuperacao de C e RMSE de DoA vs L (transmissoes), CEGO vs CONHECIDA,
% SEM interferidor (so' acoplamento real + ruido). SEM ORACLE.
% Depende no path: doa_kw_uca, estimate_C_circulant_pool,
%   estimate_C_selfcal_multishot, utils.steering_vec_uca.
% Saidas: <outDir>/multishot_frobenius_vs_L.png, multishot_rmse_doa_vs_L.png
% =========================================================================
fprintf('\n\n========== BLOCO MULTISHOT SELF-CAL (modelo real) ==========\n');

% --- ponto de operacao deste estudo ---
ms_SNR_dB    = 3;
ms_lambda    = 1e-2;        % Tikhonov no pool
ms_n_refine  = 0;          % 0 = sem refino
ms_L_list    = [1 2 4 6 8 12];
ms_n_rep     = 60;
ms_iR        = 1;          % indice de raio (range_radius); script usa 1
ms_radius_m  = range_radius(ms_iR)*lambda;
ms_theta_deg = theta_cal_deg;                 % plano XY (90)

% --- matriz de acoplamento REAL (a mesma do canal) ---
C_true_ms = Coupling_matrices_real(:,:,ms_iR);
CtN_ms  = C_true_ms / C_true_ms(1,1);
froT_ms = norm(CtN_ms, 'fro');

% --- preambulo FIXO: reusa o q_cal ja' gerado na calibracao ---
qv_ms = q_cal(:);
qv_ms = qv_ms / sqrt(mean(abs(qv_ms).^2) + eps);     % potencia 1
K_ms  = numel(qv_ms);

% --- [DIAG] sanidade ---
fprintf('[multishot] raio=%.2f lambda | ||C nao-circ||: %.3f\n', ...
        range_radius(ms_iR), norm(C_true_ms-mean(diag(C_true_ms))*eye(M),'fro'));
beta_ms = 2*pi*(0:M-1).'/M;
L_diag = 6; err_raw = zeros(1,L_diag);
for l = 1:L_diag
    phi_sig = randi([-180,179]);
    a_s = utils.steering_vec_uca(M, ms_radius_m, lambda, ms_theta_deg, phi_sig);
    Xn  = sqrt(10^(-ms_SNR_dB/10)/2)*(randn(M,K_ms)+1j*randn(M,K_ms));
    Xw  = C_true_ms*(a_s*qv_ms.') + Xn;
    [~, phk, ~] = doa_kw_uca(Xw, qv_ms.', ms_radius_m, lambda, beta_ms);
    e = abs(phk(1)-phi_sig); e = min(e,360-e); err_raw(l) = e;
end
fprintf('[multishot] DIAG RMSE KW cru = %.2f deg\n', sqrt(mean(err_raw.^2)));

% --- varredura em L ---
nL = numel(ms_L_list);
fro_ms = zeros(2,nL); doa_ms = zeros(2,nL);    % [cego; conhecida]
for il = 1:nL
  L = ms_L_list(il); acc = zeros(2,2);
  for rep = 1:ms_n_rep
    Xw = cell(1,L); phis_true = zeros(1,L);
    for l = 1:L
      phi_sig = randi([-180,179]); phis_true(l) = phi_sig;
      a_s = utils.steering_vec_uca(M, ms_radius_m, lambda, ms_theta_deg, phi_sig);
      Xn_w = sqrt(1/10^(ms_SNR_dB/10)/2)*(randn(M,K_ms)+1j*randn(M,K_ms));
      Xw{l} = C_true_ms*(a_s*qv_ms.') + Xn_w;
    end
    % CEGO
    [Cm, phim] = estimate_C_selfcal_multishot(Xw, qv_ms, M, ms_radius_m, lambda, ms_n_refine, ms_lambda, eye(M));
    em = abs(phim-phis_true); em=min(em,360-em);
    acc(1,1)=acc(1,1)+norm(Cm/Cm(1,1)-CtN_ms,'fro')/froT_ms; acc(1,2)=acc(1,2)+mean(em.^2);
    % CONHECIDA
    A_known = zeros(M,L); B_known = zeros(M,L);
    for l=1:L
      A_known(:,l) = utils.steering_vec_uca(M, ms_radius_m, lambda, ms_theta_deg, phis_true(l));
      B_known(:,l) = Xw{l}*conj(qv_ms)/(qv_ms'*qv_ms);
    end
    Ck = estimate_C_circulant_pool(B_known, A_known, M, ms_lambda, eye(M));
    acc(2,1)=acc(2,1)+norm(Ck/Ck(1,1)-CtN_ms,'fro')/froT_ms;
  end
  fro_ms(:,il)=acc(:,1)/ms_n_rep; doa_ms(:,il)=sqrt(acc(:,2)/ms_n_rep);
  fprintf('  L=%2d | Frob cego=%.3f conhecida=%.3f | RMSE_DoA=%.2f\n', ...
          L, fro_ms(1,il), fro_ms(2,il), doa_ms(1,il));
end

% --- plots ---
ms_cols = [0.13 0.47 0.71; 0.5 0.5 0.5];
fms1 = figure('Name','C-recovery vs L (modelo real)','NumberTitle','off','Position',[60 60 720 480]);
hold on; grid on;
plot(ms_L_list, fro_ms(1,:), '-s', 'Color',ms_cols(1,:),'LineWidth',1.8,'MarkerSize',6,'MarkerFaceColor',ms_cols(1,:));
plot(ms_L_list, fro_ms(2,:), '--o','Color',ms_cols(2,:),'LineWidth',1.8,'MarkerSize',6,'MarkerFaceColor',ms_cols(2,:));
xlabel('L (transmissoes bufferizadas)'); ylabel('||C^{hat}-C^{true}||_F / ||C^{true}||_F');
title(sprintf('Recuperacao de C (multishot CEGO, KW, modelo real)  |  SNR=%+d, \\lambda=%.0e', ms_SNR_dB, ms_lambda));
legend({'CEGO (direcoes via KW)','direcoes conhecidas (ref)'},'Location','northeast');
xlim([0 max(ms_L_list)+1]); ylim([0 max(0.5, max(fro_ms(:))*1.05)]);
exportgraphics(fms1, fullfile(outDir,'multishot_frobenius_vs_L.png'), 'Resolution', 200);

fms2 = figure('Name','DoA RMSE vs L (modelo real)','NumberTitle','off','Position',[80 60 720 480]);
hold on; grid on;
plot(ms_L_list, doa_ms(1,:), '-s','Color',ms_cols(1,:),'LineWidth',1.8,'MarkerSize',6,'MarkerFaceColor',ms_cols(1,:));
xlabel('L (transmissoes bufferizadas)'); ylabel('RMSE de \phi (graus)');
title(sprintf('RMSE de DoA por transmissao (KW, multishot CEGO, modelo real)  |  SNR=%+d', ms_SNR_dB));
xlim([0 max(ms_L_list)+1]);
exportgraphics(fms2, fullfile(outDir,'multishot_rmse_doa_vs_L.png'), 'Resolution', 200);
fprintf('========== FIM BLOCO MULTISHOT ==========\n\n');
