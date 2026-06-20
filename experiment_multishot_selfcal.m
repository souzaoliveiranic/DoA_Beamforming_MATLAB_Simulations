% =========================================================================
% EXPERIMENTO: Auto-calibracao CEGA multi-transmissao (KW only)
% USANDO O MODELO REAL DO SCRIPT PRINCIPAL (sem interferidor).
% -------------------------------------------------------------------------
% Matriz de acoplamento: compute_Ctx_for_R (Antenna Toolbox) + quebra de
%   circularidade de 2%, EXATAMENTE como no script principal (Coupling_matrices_real).
% Geracao de sinal: utils.simulate_fsk_data_uca (mesmo modelo FSK-2 do script).
% Preambulo FIXO reusado em todas as transmissoes (premissa do multishot).
%
% Estuda recuperacao de C e RMSE de DoA vs L (transmissoes bufferizadas),
% comparando CEGO (direcoes via KW) vs CONHECIDAS (referencia).
%
% >>> REQUER o ambiente do script principal:
%     - Antenna Toolbox (compute_Ctx_for_R usa dipole/circularArray/sparameters)
%     - utils.m no path (simulate_fsk_data_uca, steering_vec_uca)
%     - doa_kw_uca, estimate_C_circulant_pool, estimate_C_selfcal_multishot no path
%     - compute_Ctx_for_R no path (e' funcao local do seu script; copie-a para
%       um arquivo compute_Ctx_for_R.m, ou rode este bloco dentro do script).
%
% Saidas: multishot_frobenius_vs_L.png, multishot_rmse_doa_vs_L.png
% SEM ORACLE: a self-cal roda cega; C_true entra so' na medicao.
% =========================================================================
clear; close all; clc;

% ---- parametros (mesmos do script principal) ----
M        = 8;
fc       = 500e6;  c0 = 3e8;  lambda = c0/fc;
r_over_l = 0.25;   exp_radius_m = r_over_l*lambda;   % range_radius do script
fs       = 288000;
N        = 2100;                                     % amostras (=K_kw)
K_kw     = 2100;
Rs       = 9600;   sps = 30;  alpha = 0.3;  span = 8;  fd = 4.8e3;
theta_sig_deg = 90;
Z0       = 50;
C_noncirc_level = 0.02;        % 2% quebra de circulancia (como no script)
outDir   = pwd;

% ---- ponto de operacao (sem interferidor) ----
ms_SNR_dB   = 3;
ms_lambda   = 1e-2;
ms_n_refine = 0;
L_list      = [1 2 4 6 8 12];
ms_n_rep    = 60;
ms_seed_circ = 20260101;       % semente da quebra de circulancia (reprodutivel)

% =====================================================================
% MATRIZ DE ACOPLAMENTO REAL (igual ao script principal, linhas 388-415)
% =====================================================================
Ctx = compute_Ctx_for_R(fc, M, exp_radius_m, Z0);     % circulante ideal (Antenna TB)
rng(ms_seed_circ);
pert_aditiva = (randn(M) + 1j*randn(M))/sqrt(2) * C_noncirc_level * mean(abs(Ctx(:)));
C_true_ms = Ctx + pert_aditiva;                       % matriz REAL (nao-circulante)
rng('shuffle');
CtN  = C_true_ms / C_true_ms(1,1);
froT = norm(CtN, 'fro');
fprintf('Matriz de acoplamento real: ||pert nao-circ||/||Ctx|| = %.3f\n', ...
        norm(C_true_ms-Ctx,'fro')/norm(Ctx,'fro'));

% =====================================================================
% PREAMBULO FIXO: gerado UMA vez via simulate_fsk_data_uca, reusado sempre.
% (o modulador usa randi internamente, entao capturamos UMA waveform e a
%  reaplicamos em todas as transmissoes -> premissa do multishot.)
% =====================================================================
rng(20260202);
[~, q_fix, ~, ~, ~, ~] = utils.simulate_fsk_data_uca(M, exp_radius_m, lambda, ...
    0, 90, theta_sig_deg, theta_sig_deg, 100, -100, N, fs, Rs, sps, alpha, span, fd);
rng('shuffle');
qv_fixed = q_fix(:);
qv_fixed = qv_fixed / sqrt(mean(abs(qv_fixed).^2) + eps);   % potencia 1
K_kw = numel(qv_fixed);

% =====================================================================
% [DIAG] Sanidade
% =====================================================================
fprintf('\n===== [DIAG] sanidade (L=6, sem interferidor) =====\n');
L_diag = 6; Xw_d = cell(1,L_diag); phis_d = zeros(1,L_diag);
for l = 1:L_diag
    phi_sig = randi([-180,179]); phis_d(l) = phi_sig;
    a_s = utils.steering_vec_uca(M, exp_radius_m, lambda, theta_sig_deg, phi_sig);
    Xn = sqrt(10^(-ms_SNR_dB/10)/2)*(randn(M,K_kw)+1j*randn(M,K_kw));
    Xw_d{l} = C_true_ms*(a_s*qv_fixed.') + Xn;
end
beta_d = 2*pi*(0:M-1).'/M;
fprintf(' KW cru por transmissao:\n'); err_raw = zeros(1,L_diag);
for l = 1:L_diag
    [~, phi_kw, ~] = doa_kw_uca(Xw_d{l}, qv_fixed.', exp_radius_m, lambda, beta_d);
    e = abs(phi_kw(1)-phis_d(l)); e = min(e,360-e); err_raw(l) = e;
    fprintf('   l=%d: true=%4d hat=%4.0f erro=%4.0f\n', l, phis_d(l), phi_kw(1), e);
end
fprintf('   RMSE KW cru = %.2f deg\n', sqrt(mean(err_raw.^2)));
fprintf('===== fim DIAG =====\n\n');

% =====================================================================
% VARREDURA EM L
% =====================================================================
nL = numel(L_list);
fro_acc = zeros(2, nL);   % [cego; conhecida]
doa_acc = zeros(2, nL);

fprintf('[multishot] SNR=%+d lambda=%.0e raio=%.2f lambda (modelo real, sem interf)\n', ...
        ms_SNR_dB, ms_lambda, r_over_l);
for il = 1:nL
  L = L_list(il); acc = zeros(2,2);
  for rep = 1:ms_n_rep
    Xw = cell(1,L); phis_true = zeros(1,L);
    for l = 1:L
      phi_sig = randi([-180,179]); phis_true(l) = phi_sig;
      a_s = utils.steering_vec_uca(M, exp_radius_m, lambda, theta_sig_deg, phi_sig);
      Xsig_w = a_s * qv_fixed.';
      Xn_w = sqrt(1/10^(ms_SNR_dB/10)/2)*(randn(M,K_kw)+1j*randn(M,K_kw));
      Xw{l} = C_true_ms*Xsig_w + Xn_w;
    end
    % (1) CEGO
    [Cm, phim] = estimate_C_selfcal_multishot(Xw, qv_fixed, M, exp_radius_m, lambda, ms_n_refine, ms_lambda, eye(M));
    em = abs(phim-phis_true); em=min(em,360-em);
    acc(1,1)=acc(1,1)+norm(Cm/Cm(1,1)-CtN,'fro')/froT; acc(1,2)=acc(1,2)+mean(em.^2);
    % (2) CONHECIDA
    A_known = zeros(M,L); B_known = zeros(M,L);
    for l=1:L
      A_known(:,l) = utils.steering_vec_uca(M, exp_radius_m, lambda, theta_sig_deg, phis_true(l));
      B_known(:,l) = Xw{l}*conj(qv_fixed)/(qv_fixed'*qv_fixed);
    end
    Ck = estimate_C_circulant_pool(B_known, A_known, M, ms_lambda, eye(M));
    acc(2,1)=acc(2,1)+norm(Ck/Ck(1,1)-CtN,'fro')/froT;
  end
  fro_acc(:,il) = acc(:,1)/ms_n_rep;
  doa_acc(:,il) = sqrt(acc(:,2)/ms_n_rep);
  fprintf('  L=%2d | Frob cego=%.3f conhecida=%.3f | RMSE_DoA=%.2f\n', ...
          L, fro_acc(1,il), fro_acc(2,il), doa_acc(1,il));
end

% =====================================================================
% PLOTS
% =====================================================================
cols = [0.13 0.47 0.71; 0.5 0.5 0.5];
f1 = figure('Name','C-recovery vs L (modelo real)','NumberTitle','off','Position',[60 60 720 480]);
hold on; grid on;
plot(L_list, fro_acc(1,:), '-s', 'Color',cols(1,:),'LineWidth',1.8,'MarkerSize',6,'MarkerFaceColor',cols(1,:));
plot(L_list, fro_acc(2,:), '--o','Color',cols(2,:),'LineWidth',1.8,'MarkerSize',6,'MarkerFaceColor',cols(2,:));
xlabel('L (transmissoes bufferizadas)'); ylabel('||C^{hat}-C^{true}||_F / ||C^{true}||_F');
title(sprintf('Recuperacao de C (multishot CEGO, KW, modelo real)  |  SNR=%+d, \\lambda=%.0e', ms_SNR_dB, ms_lambda));
legend({'CEGO (direcoes via KW)','direcoes conhecidas (ref)'},'Location','northeast');
xlim([0 max(L_list)+1]); ylim([0 max(0.5, max(fro_acc(:))*1.05)]);
exportgraphics(f1, fullfile(outDir,'multishot_frobenius_vs_L.png'), 'Resolution', 200);
fprintf('  -> multishot_frobenius_vs_L.png\n');

f2 = figure('Name','DoA RMSE vs L (modelo real)','NumberTitle','off','Position',[80 60 720 480]);
hold on; grid on;
plot(L_list, doa_acc(1,:), '-s','Color',cols(1,:),'LineWidth',1.8,'MarkerSize',6,'MarkerFaceColor',cols(1,:));
xlabel('L (transmissoes bufferizadas)'); ylabel('RMSE de \phi (graus)');
title(sprintf('RMSE de DoA por transmissao (KW, multishot CEGO, modelo real)  |  SNR=%+d', ms_SNR_dB));
xlim([0 max(L_list)+1]);
exportgraphics(f2, fullfile(outDir,'multishot_rmse_doa_vs_L.png'), 'Resolution', 200);
fprintf('  -> multishot_rmse_doa_vs_L.png\n');
fprintf('[multishot] concluido.\n');














function Ctx = compute_Ctx_for_R(fc, M, R, Z0)
% Monta UCA e calcula Z via S-parameters, depois Ctx (pág. 29)

c = 3e8;
lambda = c/fc;

% --- Elemento ---
mp = dipole;
mp.Length  = 0.5*lambda;
mp.Width  = 0.01*lambda;

% --- UCA ---
uca = circularArray;
uca.Element = mp;
uca.NumElements = M;
uca.Radius = R;

% --- S-parameters = matriz de acoplamento ---
Sobj = sparameters(uca, fc);
S_matrix = Sobj.Parameters(:,:,1);

Z_matrix = s2z(S_matrix, Z0); % conversão Z corrijida

Zg = Z0;
Zself = diag(Z_matrix);          % Z11, Z22, ..., ZNN
denom = Zself + Zg;              % (Zjj + Zg,j)

Ctx = eye(M);

% fora-diagonais: C(i,j) = Z(i,j)/(Z(j,j)+Zg(j))
for j = 1:M
    for i = 1:M
        if i ~= j
            Ctx(i,j) = Z_matrix(i,j) / denom(j);
        end
    end
end
end


   function C = build_C_circulant_uca4(c1, c2, c3, c4)
    % UCA de 8 elementos
    first_row = [1, c1, c2, c3, c4, c3, c2, c1];
    M = 8;
    C = zeros(M,M);
    for i = 1:M
        C(i,:) = circshift(first_row, [0 i-1]);
    end
end

function J = cost_mcm_circulant4(p, a, b_hat)
    c1 = p(1) + 1j*p(2);
    c2 = p(3) + 1j*p(4);
    c3 = p(5) + 1j*p(6);
    c4 = p(7) + 1j*p(8);

    C = build_C_circulant_uca4(c1, c2, c3, c4);
    b_model = C * a;

    % remove ambiguidade de ganho escalar complexo
    alpha = (b_model' * b_hat) / (b_model' * b_model);

    err = b_hat - alpha * b_model;
    J = norm(err)^2;
end
