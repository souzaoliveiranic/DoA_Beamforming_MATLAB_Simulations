function y_shifted = apply_fractional_shift(y, shift_samples)
%APPLY_FRACTIONAL_SHIFT Aplica um deslocamento (possivelmente fracionario) a
%   um sinal complexo, usando shift de fase no dominio da frequencia.
%
%   y_shifted = apply_fractional_shift(y, shift_samples)
%
% Entradas:
%   y             : Nx1 ou 1xN, sinal complexo
%   shift_samples : escalar (pode ser fracionario)
%                   shift > 0 : atraso (sinal aparece "depois")
%                   shift < 0 : adianto
%
% Saida:
%   y_shifted     : sinal deslocado, mesmo tamanho de y
%
% Modela erro de sincronismo temporal de fracao de amostra, tipico em
% receptores SDR antes do passo de timing recovery. Usa propriedade da DFT:
% multiplicar por e^{-j 2*pi*k*shift/N} no dominio da frequencia equivale a
% deslocar o sinal por "shift" amostras no dominio do tempo (com extensao
% periodica).

    if shift_samples == 0
        y_shifted = y;
        return;
    end

    is_row = isrow(y);
    y = y(:);
    N = numel(y);

    % Frequencias normalizadas centradas em zero (para shift fracionario correto)
    if mod(N, 2) == 0
        k = [(0:N/2-1), (-N/2:-1)].';
    else
        k = [(0:(N-1)/2), -(N-1)/2:-1].';
    end

    Y = fft(y);
    phase_shift = exp(-1j * 2*pi * k * shift_samples / N);
    y_shifted = ifft(Y .* phase_shift);

    if is_row
        y_shifted = y_shifted.';
    end
end
