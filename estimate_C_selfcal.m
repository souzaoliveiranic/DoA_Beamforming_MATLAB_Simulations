function [C_hat, c_hat, alpha_hat, phi_hat_deg, n_iter] = ...
    estimate_C_selfcal(X, q, M, radius, lambda, phi_grid_deg, max_iter, tol_phi, tol_C)
% ESTIMATE_C_SELFCAL  Self-calibration: estima conjuntamente DoA e matriz
%   de acoplamento usando APENAS o conhecimento da forma de onda q (KW).
%   Nao requer direcao de calibracao conhecida.
%
%   [C_hat, c_hat, alpha_hat, phi_hat_deg, n_iter] = ...
%       estimate_C_selfcal(X, q, M, radius, lambda, phi_grid_deg, ...
%                          max_iter, tol_phi, tol_C)
%
% Entradas:
%   X            : MxN, dados recebidos (1 fonte util, possivel interferente)
%   q            : Nx1 ou 1xN, forma de onda do sinal util (conhecida)
%   M            : numero de elementos
%   radius       : raio do UCA (em metros)
%   lambda       : comprimento de onda (em metros)
%   phi_grid_deg : grid de busca para o KW-DoA  (padrao -180:0.5:180)
%   max_iter     : max. iteracoes alternantes  (padrao 15)
%   tol_phi      : tolerancia em |phi_t - phi_{t-1}| (graus)  (padrao 0.01)
%   tol_C        : tolerancia em ||C - C_prev||_F / ||C_prev||_F  (padrao 1e-8)
%
% Algoritmo (alternante):
%   Init: C = I_M, phi inicial = KW-DoA(X, q)
%   Repeat:
%     1) D = inv(C);  Y = D*X
%     2) phi <- argmax_phi |a(phi)' * (Y q* / q'q)|^2 / ||a(phi)||^2
%     3) C <- estimate_C_circulant_uca(X q* / q'q, a(phi), M)
%   Until convergencia em phi e C.
%
% Vantagens:
%   - Nao precisa fase de calibracao com direcao conhecida.
%   - Usa o proprio sinal de operacao (toda transmissao "calibra").
%
% Limitacoes:
%   - Convergencia local: depende do init via KW-DoA com C=I (ja normalmente
%     bom porque acoplamento desloca o pico mas raramente o desvia muito).
%   - 1 fonte util. Para multiplas fontes, usa-se versao MUSIC-like
%     (nao implementada aqui).

    if nargin < 6 || isempty(phi_grid_deg), phi_grid_deg = -180:0.5:180; end
    if nargin < 7 || isempty(max_iter),     max_iter = 15;               end
    if nargin < 8 || isempty(tol_phi),      tol_phi  = 0.01;             end
    if nargin < 9 || isempty(tol_C),        tol_C    = 1e-8;             end

    q = q(:);
    qHq = q' * q;

    % b_hat sem compensacao (so calculado uma vez; muda apenas D*X no DoA)
    b_hat_orig = X * conj(q) / qHq;

    C_hat = eye(M);
    phi_prev = NaN;
    C_prev   = C_hat;
    n_iter   = 0;

    for it = 1:max_iter
        % --- 1. Aplica D = inv(C) e calcula b_hat compensado p/ DoA ---
        D = inv(C_hat);
        b_hat_comp = D * b_hat_orig;     % equivalente a (D*X) q*/(q'q)

        % --- 2. KW-DoA grid search ---
        best_score = -inf;
        best_phi   = NaN;
        for ip = 1:numel(phi_grid_deg)
            a_p = utils.steering_vec_uca(M, radius, lambda, 90, phi_grid_deg(ip));
            score = abs(a_p' * b_hat_comp)^2 / real(a_p' * a_p);
            if score > best_score
                best_score = score;
                best_phi   = phi_grid_deg(ip);
            end
        end
        phi_hat_deg = best_phi;

        % --- 3. Re-estima C com b_hat_orig (NAO compensado!) e a(phi_hat) ---
        a_hat = utils.steering_vec_uca(M, radius, lambda, 90, phi_hat_deg);
        [C_new, c_hat, alpha_hat, ~] = estimate_C_circulant_uca(b_hat_orig, a_hat, M);

        % --- Criterio de parada ---
        d_phi = abs(phi_hat_deg - phi_prev);
        d_C   = norm(C_new - C_prev, 'fro') / max(norm(C_prev,'fro'), eps);
        n_iter = it;

        C_hat    = C_new;
        phi_prev = phi_hat_deg;
        C_prev   = C_new;

        if (it >= 2) && (d_phi < tol_phi) && (d_C < tol_C)
            break;
        end
    end
end
