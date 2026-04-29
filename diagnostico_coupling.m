%% DIAGNOSTICO: Verificar passo a passo o comportamento da matriz de acoplamento
% Rode isso e me mostre os plots e a saida do console
clear; clc; close all;

%% Parametros
M  = 8;
fc = 500e6;
c  = 3e8;
lambda = c / fc;
Z0 = 50;

R_list_lambda = [0.1 0.2 0.25 0.375 0.5];
R_list_m      = R_list_lambda * lambda;

fprintf('=== DIAGNOSTICO DA MATRIZ DE ACOPLAMENTO ===\n');
fprintf('f = %.0f MHz | lambda = %.2f m | M = %d\n\n', fc/1e6, lambda, M);

%% 1) Verificar Z_matrix do Antenna Toolbox para cada raio
figure('Position', [50 50 1400 900], 'Color', 'w');

for k = 1:numel(R_list_m)
    R = R_list_m(k);
    
    % --- Monta UCA ---
    mp = dipole;
    mp.Length = 0.5*lambda;
    mp.Width  = 0.01*lambda;
    uca = circularArray;
    uca.Element = mp;
    uca.NumElements = M;
    uca.Radius = R;
    
    % --- S -> Z ---
    Sobj = sparameters(uca, fc);
    S_matrix = Sobj.Parameters(:,:,1);
    Z_matrix = s2z(S_matrix, Z0);
    
    % --- Diagnostico basico ---
    Z11 = Z_matrix(1,1);
    Z12 = Z_matrix(1,2);
    
    % Distancia entre adjacentes
    d_adj = 2*R*sin(pi/M);
    
    fprintf('--- R = %.3f lambda (%.4f m) ---\n', R_list_lambda(k), R);
    fprintf('  d_adjacente = %.4f lambda (%.4f m)\n', d_adj/lambda, d_adj);
    fprintf('  Z_11 = %.2f + j%.2f  (|Z_11| = %.2f)\n', real(Z11), imag(Z11), abs(Z11));
    fprintf('  Z_12 = %.2f + j%.2f  (|Z_12| = %.2f)\n', real(Z12), imag(Z12), abs(Z12));
    fprintf('  |Z_12/Z_11| = %.4f\n', abs(Z12)/abs(Z11));
    fprintf('  |Z_12|/Z0   = %.4f\n', abs(Z12)/Z0);
    
    % --- Tres formulas para C ---
    Zself = diag(Z_matrix);
    
    % Formula A: Codigo original (pag 29 TX, sem sinal negativo)
    C_A = eye(M);
    for i = 1:M
        for j = 1:M
            if i ~= j
                C_A(i,j) = Z_matrix(i,j) / (Zself(j) + Z0);
            end
        end
    end
    
    % Formula B: Com sinal negativo (pag 43 RX)
    C_B = eye(M);
    for i = 1:M
        for j = 1:M
            if i ~= j
                C_B(i,j) = -Z_matrix(i,j) / (Zself(j) + Z0);
            end
        end
    end
    
    invC_A = inv(C_A);
    invC_B = C_A; %inv(C_B);
    
    fprintf('  cond(C_A) = %.2f | cond(C_B) = %.2f\n', cond(C_A), cond(C_B));
    fprintf('  |C_A(1,2)| = %.4f | |C_B(1,2)| = %.4f\n', abs(C_A(1,2)), abs(C_B(1,2)));
    fprintf('  |invC_A(1,2)| = %.4f | |invC_B(1,2)| = %.4f\n', abs(invC_A(1,2)), abs(invC_B(1,2)));
    fprintf('\n');
    
    % --- Plot Linha 1: |Z(1,:)| para cada raio ---
    subplot(3, numel(R_list_m), k);
    bar(abs(Z_matrix(1,:)));
    title(sprintf('|Z_{1,j}| | R=%.2f\\lambda', R_list_lambda(k)));
    xlabel('j'); ylabel('Ohm');
    
    % --- Plot Linha 2: |C(:,1)| DIRETO (sem inverter) ---
    subplot(3, numel(R_list_m), numel(R_list_m) + k);
    plot(1:M, abs(C_A(:,1)), 'bo-', 'LineWidth', 1.5); hold on;
    plot(1:M, abs(C_B(:,1)), 'rs-', 'LineWidth', 1.5); hold off;
    title(sprintf('|C(:,1)| direto | R=%.2f\\lambda', R_list_lambda(k)));
    xlabel('Elemento i'); ylabel('|C_{i,1}|');
    if k == 1
        legend('A: +Z/(Z_{jj}+Z_0)', 'B: -Z/(Z_{jj}+Z_0)', 'Location', 'best');
    end
    
    % --- Plot Linha 3: |inv(C)(:,1)| ---
    e1 = zeros(M,1); e1(1) = 1;
    subplot(3, numel(R_list_m), 2*numel(R_list_m) + k);
    plot(1:M, abs(invC_A * e1), 'bo-', 'LineWidth', 1.5); hold on;
    plot(1:M, abs(invC_B * e1), 'rs-', 'LineWidth', 1.5); hold off;
    title(sprintf('|inv(C) e_1| | R=%.2f\\lambda', R_list_lambda(k)));
    xlabel('Elemento i'); ylabel('|inv(C)_{i,1}|');
end

sgtitle('Diagnostico: Z, C direto, e inv(C) para cada raio', 'FontSize', 14);

%% 2) Resumo: |C(1,2)| e |inv(C)(1,2)| vs raio
figure('Position', [100 100 900 400], 'Color', 'w');

C12_B    = zeros(1, numel(R_list_m));
invC12_B = zeros(1, numel(R_list_m));
Z12_all  = zeros(1, numel(R_list_m));
cond_B   = zeros(1, numel(R_list_m));

for k = 1:numel(R_list_m)
    R = R_list_m(k);
    mp = dipole; mp.Length = 0.5*lambda; mp.Width = 0.01*lambda;
    uca = circularArray;
    uca.Element = mp; uca.NumElements = M; uca.Radius = R;
    
    Sobj = sparameters(uca, fc);
    S_matrix = Sobj.Parameters(:,:,1);
    Z_matrix = s2z(S_matrix, Z0);
    Zself = diag(Z_matrix);
    
    Z12_all(k) = abs(Z_matrix(1,2));
    
    C_B = eye(M);
    for i = 1:M
        for j = 1:M
            if i ~= j
                C_B(i,j) = -Z_matrix(i,j) / (Zself(j) + Z0);
            end
        end
    end
    
    C12_B(k) = abs(C_B(1,2));
    tmpB = inv(C_B);
    invC12_B(k) = abs(tmpB(1,2));
    cond_B(k) = cond(C_B);
end

subplot(1,3,1);
plot(R_list_lambda, Z12_all, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'k');
xlabel('R (\lambda)'); ylabel('|Z_{12}| (Ohm)');
title('|Z_{12}| adjacente vs Raio'); grid on;

subplot(1,3,2);
plot(R_list_lambda, C12_B, 'rs-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xlabel('R (\lambda)'); ylabel('|C(1,2)|');
title('|C(1,2)| direto vs Raio'); grid on;

subplot(1,3,3);
plot(R_list_lambda, invC12_B, 'bs-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('R (\lambda)'); ylabel('|inv(C)(1,2)|');
title('|inv(C)(1,2)| vs Raio'); grid on;

sgtitle('Resumo: como acoplamento varia com o raio (Formula B)', 'FontSize', 14);

fprintf('=== FIM DO DIAGNOSTICO ===\n');