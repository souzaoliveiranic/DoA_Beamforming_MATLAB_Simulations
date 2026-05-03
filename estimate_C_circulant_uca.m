function [C_hat, c_hat, alpha_hat, residual] = estimate_C_circulant_uca(b_hat, a, M)
% ESTIMATE_C_CIRCULANT_UCA  Estimação LS da matriz de acoplamento de um UCA,
%   explorando a estrutura circulante simétrica.
%
%   [C_hat, c_hat, alpha_hat, residual] = estimate_C_circulant_uca(b_hat, a, M)
%
% Modela:
%       b_hat = alpha * C * a + ruido
% onde C é circulante com primeira linha:
%   M par   :  [1, c_1, c_2, ..., c_{K-1}, c_K, c_{K-1}, ..., c_2, c_1]   (K = M/2)
%   M ímpar :  [1, c_1, c_2, ..., c_K,            c_K, ..., c_2, c_1]    (K = (M-1)/2)
%
% Resolve por mínimos quadrados lineares (1 linha de código), explorando que
%       b_hat[i] = alpha * a[i] + sum_{k=1..K} gamma_k * S_k[i]
% com gamma_k = alpha * c_k e
%       S_k[i] = a[i+k] + a[i-k]      se k < M/2
%       S_k[i] = a[i+K]               se k = M/2  (somente quando M é par)
% (índices cíclicos, módulo M).
%
% Total: 1 + K incógnitas complexas (alpha + K coef. de acoplamento únicos),
%        M equações. Com M >= 1 + K (o que vale sempre que M >= 2),
%        sistema sobredeterminado -> LS linear.
%
% Entradas:
%   b_hat : Mx1, assinatura espacial estimada com acoplamento
%   a     : Mx1, vetor de direção ideal (sem acoplamento) na direção de calibração
%   M     : escalar, número de elementos do UCA
%
% Saídas:
%   C_hat     : MxM, matriz de acoplamento estimada (estrutura circulante)
%   c_hat     : Kx1, coeficientes de acoplamento únicos (c_1,...,c_K)
%   alpha_hat : escalar complexo, ganho global (resolve ambiguidade de escala)
%   residual  : Mx1, resíduo do ajuste LS (b_hat - alpha_hat * C_hat * a)
%
% Referência:
%   S. Khan, H. Sajjad, M. K. Ozdemir e E. Arvas,
%   "Mutual Coupling Compensation in Receiving Antenna Arrays,"
%   2020 ACES Symposium.

    b_hat = b_hat(:);
    a     = a(:);
    if numel(b_hat) ~= M || numel(a) ~= M
        error('estimate_C_circulant_uca:dim', ...
              'b_hat e a devem ser vetores Mx1.');
    end

    K = floor(M/2);
    is_even = (mod(M, 2) == 0);

    % ----- Constrói matriz LS  M_ls  de dimensão M x (1+K) ---------------
    % Coluna 1   : a[i]            (multiplica alpha)
    % Coluna 1+k : S_k[i]          (multiplica gamma_k = alpha*c_k)
    M_ls = zeros(M, 1 + K);
    for i = 1:M
        M_ls(i, 1) = a(i);
        for k = 1:K
            i_plus  = mod(i-1+k, M) + 1;
            i_minus = mod(i-1-k, M) + 1;
            if is_even && (k == K)
                % distância máxima M/2: +k e -k coincidem -> contrib única
                M_ls(i, 1+k) = a(i_plus);
            else
                M_ls(i, 1+k) = a(i_plus) + a(i_minus);
            end
        end
    end

    % ----- Resolve LS:  theta = argmin || M_ls * theta - b_hat ||^2 -------
    theta     = M_ls \ b_hat;
    alpha_hat = theta(1);
    gamma_hat = theta(2:end);

    if abs(alpha_hat) < 1e-12
        warning('estimate_C_circulant_uca:smallAlpha', ...
            'alpha_hat ~ 0 (%.3e). LS pode estar mal-condicionado.', ...
            abs(alpha_hat));
    end

    % Coeficientes de acoplamento normalizados
    c_hat = gamma_hat / alpha_hat;

    % ----- Reconstrói C_hat circulante ------------------------------------
    if is_even
        first_row = [1, c_hat(1:K-1).', c_hat(K), flip(c_hat(1:K-1).')];
    else
        first_row = [1, c_hat(1:K).', flip(c_hat(1:K).')];
    end

    C_hat = zeros(M, M);
    for i = 1:M
        C_hat(i, :) = circshift(first_row, [0, i-1]);
    end

    residual = b_hat - alpha_hat * C_hat * a;
end
