% test_estimacao_methods.m
% Validacao numerica e comparacao dos 4 estimadores de C circulante:
%   - LS 1 direcao (referencia, estimate_C_circulant_uca)
%   - Davies modal (estimate_C_davies)
%   - Multi-direcao alternante (estimate_C_multidir)
%   - Self-cal alternante (estimate_C_selfcal)

clear; close all; clc;

%% Parametros
M  = 8;
fc = 2.4e9;
c  = 3e8;
lambda = c/fc;
r  = 0.5*lambda;
Z0 = 50;

% C_true a partir das Z_t do usuario
Zt_vals = [ ...
   -0.08 - 1j*11.77;
    7.19 + 1j*2.12;
   -0.06 + 1j*5.96;
   -2.50 + 1j*5.00];
c_true = -Zt_vals / Z0;
first_row_true = [1, c_true(1), c_true(2), c_true(3), ...
                     c_true(4), c_true(3), c_true(2), c_true(1)];
C_true = zeros(M, M);
for i = 1:M
    C_true(i, :) = circshift(first_row_true, [0, i-1]);
end
norm_Ctrue = norm(C_true, 'fro');

%% Helpers
beta = 2*pi*(0:M-1).'/M;
k0   = 2*pi/lambda;
steer = @(phi_deg) ...
    exp(1j*k0*r * cos(deg2rad(phi_deg) - beta)) ...
    / norm(exp(1j*k0*r * cos(deg2rad(phi_deg) - beta)));

err_F = @(C) norm(C - C_true, 'fro') / norm_Ctrue;

%% --- Caso A: SEM RUIDO ---
fprintf('==========================================================\n');
fprintf('CASO A: sem ruido\n');
fprintf('==========================================================\n');

alpha_true = 2.5*exp(1j*0.7);
phi_cal = 0;
a_cal = steer(phi_cal);
b_clean = alpha_true * (C_true * a_cal);

% LS 1-dir
[C_ls, ~, alpha_ls, ~] = estimate_C_circulant_uca(b_clean, a_cal, M);

% Davies
[C_dav, lam_dav, alpha_dav, ~] = estimate_C_davies(b_clean, a_cal, M);

% Multi-dir P=4
phis_cal = [0 30 60 90];
P = numel(phis_cal);
B_multi = zeros(M, P); A_multi = zeros(M, P);
alphas_true_p = alpha_true * exp(1j*0.3*(0:P-1)');
for p = 1:P
    a_p = steer(phis_cal(p));
    A_multi(:, p) = a_p;
    B_multi(:, p) = alphas_true_p(p) * (C_true * a_p);
end
[C_md, ~, alphas_md, ~, n_iter_md] = estimate_C_multidir(B_multi, A_multi, M);

fprintf('  LS 1-dir         : err Frob = %.3e\n', err_F(C_ls));
fprintf('  Davies modal     : err Frob = %.3e\n', err_F(C_dav));
fprintf('  Multi-dir P=%d (%d it.) : err Frob = %.3e\n', P, n_iter_md, err_F(C_md));
fprintf('  -> Esperado: todos ~1e-15 (precisao de maquina)\n\n');

%% --- Caso B: COM RUIDO ---
fprintf('==========================================================\n');
fprintf('CASO B: com ruido vetorial sigma=1e-3, 200 trials\n');
fprintf('==========================================================\n');
sigma = 1e-3;
nTr = 200;
errs_ls  = zeros(nTr,1);
errs_dav = zeros(nTr,1);
errs_md  = zeros(nTr,1);

for t = 1:nTr
    n1 = sigma * (randn(M,1) + 1j*randn(M,1)) / sqrt(2);
    b_n = b_clean + n1;
    [Cl,~,~,~] = estimate_C_circulant_uca(b_n, a_cal, M);
    [Cd,~,~,~] = estimate_C_davies(b_n, a_cal, M);
    errs_ls(t)  = err_F(Cl);
    errs_dav(t) = err_F(Cd);

    % Multi-dir
    Bn = zeros(M,P); An = zeros(M,P);
    for p = 1:P
        a_p = steer(phis_cal(p));
        n_p = sigma * (randn(M,1) + 1j*randn(M,1)) / sqrt(2);
        An(:,p) = a_p;
        Bn(:,p) = alphas_true_p(p) * (C_true * a_p) + n_p;
    end
    [Cmd,~,~,~,~] = estimate_C_multidir(Bn, An, M);
    errs_md(t) = err_F(Cmd);
end
fprintf('  LS 1-dir         : mediana %.3e | max %.3e\n', median(errs_ls),  max(errs_ls));
fprintf('  Davies modal     : mediana %.3e | max %.3e\n', median(errs_dav), max(errs_dav));
fprintf('  Multi-dir P=%d    : mediana %.3e | max %.3e\n', P, median(errs_md), max(errs_md));
fprintf('  -> Esperado: Multi-dir ~sqrt(P) melhor que 1-dir\n\n');

%% --- Caso C: CONDICIONAMENTO Davies vs phi_cal ---
fprintf('==========================================================\n');
fprintf('CASO C: Davies modal - min|F a| vs phi_cal\n');
fprintf('==========================================================\n');
fprintf('  phi_cal | min|F a| | cond efetivo (max/min)\n');
fprintf('  --------+----------+-----------------------\n');
for phi_test = [0 5 10 15 20 22.5 30 45]
    a_t = steer(phi_test);
    Fa = fft(a_t);
    fprintf('  %6.1f  | %.3e | %.3e\n', phi_test, min(abs(Fa)), max(abs(Fa))/min(abs(Fa)));
end
fprintf('  Atencao: phi_cal=22.5 deg fica entre antenas -> modo cego.\n');
fprintf('  Use phi_cal=0 deg (cond ~1.78, otimo).\n');
