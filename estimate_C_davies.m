function [C_hat, lambda_hat, alpha_hat, residual] = estimate_C_davies(b_hat, a, M)
% ESTIMATE_C_DAVIES  Estimacao da matriz de acoplamento via diagonalizacao DFT.
%
%   [C_hat, lambda_hat, alpha_hat, residual] = estimate_C_davies(b_hat, a, M)
%
% Modela:    b_hat = alpha * C * a + ruido,  C circulante MxM.
% Como toda matriz circulante e' diagonalizada pela DFT,
%      F * C * F^H = diag(lambda),     com F = dftmtx(M),
% temos no dominio modal:
%      (F * b_hat)_m = alpha * lambda_m * (F * a)_m
% Logo cada bin m da uma equacao escalar independente:
%      mu_m = (F b)_m / (F a)_m   =   alpha * lambda_m
% E como o nosso modelo impoe c_0 = 1 (primeiro elemento da linha circulante),
% temos sum(lambda)/M = c_0 = 1, logo  alpha = mean(mu)  e lambda_m = mu_m/alpha.
% A primeira linha de C reconstroi-se por  ifft(lambda).
%
% Vantagens vs LS direto:
%   - 5 linhas de codigo, sem matriz LS, sem inversao.
%   - Cada modo e' independente: ve-se exatamente quais sao
%     mal-condicionados (|F a|_m proximo de zero).
%
% Cuidado: bins onde |(F a)_m| e' pequeno ficam dominados por ruido.
%   Para UCA-8 com phi_cal = 22.5 graus (entre antenas) algum bin se anula
%   exatamente -> sistema singular (mesmo problema do LS).
%   Com phi_cal = 0 graus, todos os bins tem magnitude ~ 0.5..0.8 (otimo).

    b_hat = b_hat(:);
    a     = a(:);

    if numel(b_hat) ~= M || numel(a) ~= M
        error('estimate_C_davies:dim', 'b_hat e a devem ser vetores Mx1.');
    end

    Fb = fft(b_hat);
    Fa = fft(a);

    % Aviso se algum modo for cego
    min_mag = min(abs(Fa));
    if min_mag < 1e-10
        warning('estimate_C_davies:blindMode', ...
            'Modo cego detectado: min|F a| = %.3e. Resultado mal-condicionado.', ...
            min_mag);
    end

    mu         = Fb ./ Fa;             % alpha * lambda_m
    alpha_hat  = mean(mu);             % pois c_0 = 1 => mean(lambda) = 1
    lambda_hat = mu / alpha_hat;
    first_row  = ifft(lambda_hat);     % primeira linha de C

    C_hat = zeros(M, M);
    for i = 1:M
        C_hat(i, :) = circshift(first_row.', [0, i-1]);
    end

    residual = b_hat - alpha_hat * (C_hat * a);
end
