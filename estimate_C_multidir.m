function [C_hat, c_hat, alpha_hat, residual, n_iter] = ...
    estimate_C_multidir(B_hat, A, M, max_iter, tol)
% ESTIMATE_C_MULTIDIR  Estimacao de C circulante a partir de P direcoes
%   conhecidas, alternancia entre alpha_p (1 escalar/medida) e c (vetor global).
%
%   [C_hat, c_hat, alpha_hat, residual, n_iter] = ...
%       estimate_C_multidir(B_hat, A, M, max_iter, tol)
%
% Entradas:
%   B_hat   : MxP, P assinaturas espaciais (uma por direcao de calibracao)
%   A       : MxP, P vetores de direcao ideais
%   M       : numero de elementos do UCA
%   max_iter: max. iteracoes (padrao 30)
%   tol     : tolerancia em ||C - C_prev||_F / ||C_prev||_F  (padrao 1e-10)
%
% Modelo:   b_p = alpha_p * C * a_p + n_p,  p = 1..P
%           C circulante, alpha_p especifico de cada medida
%           (waveforms diferentes ou medidas em momentos distintos).
%
% Algoritmo:
%   Init: usa primeira direcao com estimate_C_circulant_uca (LS one-shot).
%   Loop alternante:
%     (1) Dado C, estima cada alpha_p (LS escalar):
%             alpha_p = (C a_p)' b_p / ||C a_p||^2
%     (2) Dado todos alpha_p, estima c global empilhando as P*M equacoes:
%             b_p[i] - alpha_p * a_p[i]  =  alpha_p * sum_k c_k * S_pk[i]
%   Ganho de variancia: para P direcoes "bem-condicionadas",
%   var(c_hat) ~ var(c_hat_1dir) / P.

    if nargin < 4 || isempty(max_iter), max_iter = 30; end
    if nargin < 5 || isempty(tol),      tol      = 1e-10; end

    [Mb, P] = size(B_hat);
    [Ma, Pa] = size(A);
    if Mb ~= M || Ma ~= M || P ~= Pa
        error('estimate_C_multidir:dim', ...
              'B_hat e A devem ser MxP com mesmo P.');
    end

    K = floor(M/2);
    is_even = (mod(M, 2) == 0);

    % --- Init: LS one-shot na 1a direcao ---
    [C_hat, c_hat, alpha_1, ~] = estimate_C_circulant_uca(B_hat(:,1), A(:,1), M);
    alpha_hat = zeros(P, 1);
    alpha_hat(1) = alpha_1;

    C_prev = C_hat;
    n_iter = 0;
    for it = 1:max_iter
        % --- Passo 1: estima alpha_p para cada direcao ---
        for p = 1:P
            Ca = C_hat * A(:, p);
            alpha_hat(p) = (Ca' * B_hat(:, p)) / (Ca' * Ca);
        end

        % --- Passo 2: estima c global a partir de todas as P*M eqs ---
        % b_p[i] - alpha_p * a_p[i] = alpha_p * sum_k c_k * S_pk[i]
        M_big = zeros(M*P, K);
        rhs   = zeros(M*P, 1);
        for p = 1:P
            a_p = A(:, p); b_p = B_hat(:, p); ap = alpha_hat(p);
            for i = 1:M
                row_idx = (p-1)*M + i;
                for k = 1:K
                    ip = mod(i-1+k, M) + 1;
                    im = mod(i-1-k, M) + 1;
                    if is_even && (k == K)
                        M_big(row_idx, k) = ap * a_p(ip);
                    else
                        M_big(row_idx, k) = ap * (a_p(ip) + a_p(im));
                    end
                end
                rhs(row_idx) = b_p(i) - ap * a_p(i);
            end
        end
        c_hat = M_big \ rhs;

        % --- Reconstroi C ---
        if is_even
            first_row = [1, c_hat(1:K-1).', c_hat(K), flip(c_hat(1:K-1).')];
        else
            first_row = [1, c_hat(1:K).', flip(c_hat(1:K).')];
        end
        C_new = zeros(M, M);
        for i = 1:M
            C_new(i, :) = circshift(first_row, [0, i-1]);
        end

        delta = norm(C_new - C_prev, 'fro') / max(norm(C_prev,'fro'), eps);
        C_hat = C_new;
        C_prev = C_new;
        n_iter = it;
        if delta < tol
            break;
        end
    end

    % Residuo total
    residual = zeros(M*P, 1);
    for p = 1:P
        residual((p-1)*M + (1:M)) = B_hat(:,p) - alpha_hat(p) * (C_hat * A(:,p));
    end
end
