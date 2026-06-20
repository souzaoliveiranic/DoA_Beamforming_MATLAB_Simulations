function [C, c_unique] = estimate_C_circulant_pool(B, A, M, lambda, C_prior)
% ESTIMATE_C_CIRCULANT_POOL  Estima UMA matriz de acoplamento circulante a
%   partir de L assinaturas, cada uma de uma direcao DIFERENTE (steering em A).
%   Formulacao estilo Zhang-Lu-Hui (2006): LS conjunto SEM ganho por
%   observacao -- a assinatura b_l e o steering a_l estao na MESMA escala
%   (forma de onda conhecida e normalizada), entao NAO ha alpha_l a estimar.
%   Introduzir alpha_l (versao anterior) so' adicionava graus de liberdade
%   espurios que absorviam ruido e desestabilizavam a estimacao.
%
%   Modelo:  B(:,l) ~ C * A(:,l),  l=1..L   (C circulante, igual p/ todo l)
%   Beamspace (C diagonal):  Bh(p,l) = ctil_p * Ah(p,l),  ctil = fft(c)
%   LS por modo (Tikhonov p/ os modos cegos J_p(kr)~0):
%
%     ctil_p = [ sum_l conj(Ah(p,l)) Bh(p,l) + lambda*ctil_prior_p ]
%              / [ sum_l |Ah(p,l)|^2          + lambda ]
%
%   - modos visiveis (|Ah|^2 grande): a diversidade de fase e^{jp phi_l}
%     das L direcoes acumula SNR e cancela interferente movel.
%   - modos cegos (|Ah|^2 ~ 0): vao para o prior.
%
% Entradas:
%   B       : M x L assinaturas (b_l = X_l*conj(q)/(q'q)). Normalizadas aqui.
%   A       : M x L steering vectors a(phi_l) (direcoes conhecidas/estimadas).
%   M       : numero de elementos.
%   lambda  : (default 0) regularizacao de Tikhonov por modo.
%   C_prior : (default I) prior MxM (modos cegos puxados p/ ele).
%
% Saidas: C (MxM circulante, diag=1), c_unique (convencao do codigo).

    if nargin < 4 || isempty(lambda),  lambda  = 0;        end
    if nargin < 5 || isempty(C_prior), C_prior = eye(M);   end

    L = size(B, 2);
    % normaliza cada assinatura p/ norma unitaria E cada steering tambem,
    % garantindo que b_l e a_l fiquem na MESMA escala (premissa Zhang).
    for l = 1:L
        nb = norm(B(:,l));  if nb > eps, B(:,l) = B(:,l)/nb; end
        na = norm(A(:,l));  if na > eps, A(:,l) = A(:,l)/na; end
    end
    Ah = fft(A);                       % M x L  (modos de cada direcao)
    Bh = fft(B);                       % M x L
    ctil_prior = fft(C_prior(:,1));    % p/ prior=I -> vetor de uns

    % LS por modo, fechado (sem alternancia, sem alpha)
    num  = sum(conj(Ah) .* Bh, 2) + lambda * ctil_prior;
    den  = sum(abs(Ah).^2, 2)     + lambda;
    ctil = num ./ den;
    c = ifft(ctil);
    if abs(c(1)) < eps, c(1) = eps; end
    c = c / c(1);                                 % normaliza diag=1
    C = zeros(M);
    for col = 1:M, C(:,col) = circshift(c, col-1); end
    c_unique = [C(1, 2:M/2).' ; C(1, M/2+1)];
end
