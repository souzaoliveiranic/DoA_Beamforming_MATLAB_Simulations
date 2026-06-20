function [C_hat, phis_hat_deg, history] = estimate_C_selfcal_multishot( ...
    X_cell, q, M, radius, lambda, n_refine, selfcal_lambda, selfcal_prior, C_true_for_log, use_safeguard)
% ESTIMATE_C_SELFCAL_MULTISHOT  Auto-calibracao CEGA multi-transmissao (KW).
%   On-the-fly: NAO conhece as direcoes. As direcoes vem do KW aplicado aos
%   dados CRUS (C=I) -- o KW e' robusto ao acoplamento, entao NAO precisa de
%   compensacao para estimar direcao (compensar com C estimado, alias, PIORA
%   o KW; ver estudo single-shot). A matriz C e' obtida juntando as L
%   assinaturas via pool por modo (estimate_C_circulant_pool). A diversidade
%   natural de direcao (fonte em movimento) torna a estimacao de C bem-posta
%   nos modos visiveis, SEM conhecer as direcoes.
%
%   IMPORTANTE: as direcoes reportadas (phis_hat_deg) sao as do KW CRU, que
%   sao as melhores (a compensacao nao as melhora). O C serve para compensar
%   o array no processamento DOWNSTREAM (ex.: beamforming), nao para o DoA.
%
% Entradas:
%   X_cell        : cell {1 x L}, cada um MxN (janela do preambulo da transmissao l).
%   q             : Nx1, forma de onda conhecida (mesma p/ todas as transmissoes).
%   M, radius, lambda : geometria (lambda = comprimento de onda).
%   n_refine      : (default 0) passos EXTRA de refino (compensa->KW->re-pool).
%                   Default 0 = sem refino (mais robusto). Use >0 com cautela:
%                   pode degradar se o C intermediario for ruim (L pequeno).
%   selfcal_lambda: (default 1e-2) Tikhonov no pool (prende modos cegos ao prior).
%   selfcal_prior : (default I) prior MxM.
%   C_true_for_log: (opcional) SO p/ log de erro Frobenius. Sem oracle.
%   use_safeguard : (default true) se true, cada passo de refino so' e' aceito
%                   se o residuo observavel diminuir (anti-degradacao); se
%                   false, aplica todos os n_refine passos cegamente (para
%                   estudar a importancia do criterio de parada).
%
% Saidas:
%   C_hat        : MxM circulante estimada (diag=1).
%   phis_hat_deg : 1xL direcoes estimadas (KW cru, robusto).
%   history      : .errF (se C_true dado; log), .phis_raw, .c_col,
%                  .res_per_iter, .stopped_at_iter

    if nargin < 6  || isempty(n_refine),       n_refine = 0;           end
    if nargin < 7  || isempty(selfcal_lambda), selfcal_lambda = 1e-2;  end
    if nargin < 8  || isempty(selfcal_prior),  selfcal_prior = eye(M); end
    if nargin < 9,                             C_true_for_log = [];    end
    if nargin < 10 || isempty(use_safeguard),  use_safeguard = true;   end

    q = q(:);  qHq = q' * q;
    L = numel(X_cell);
    beta_uca = 2*pi*(0:M-1).'/M;

    % [DIAG] Confirma no terminal que ESTA versao (KW cru + pool, sem feedback)
    % esta carregada. Se voce NAO ver esta linha, o MATLAB tem versao antiga em
    % cache -> rode 'clear functions'. Imprime uma vez por sessao.
    persistent ms_versao_avisada
    if isempty(ms_versao_avisada)
        fprintf(['[multishot] versao SEM-FEEDBACK carregada ' ...
                 '(KW cru p/ direcoes + pool; n_refine=%d).\n'], n_refine);
        ms_versao_avisada = true;
    end

    % Assinaturas originais (uma por transmissao), fixas
    B_orig = zeros(M, L);
    for l = 1:L, B_orig(:,l) = X_cell{l} * conj(q) / qHq; end

    % --- 1. Direcoes via KW nos dados CRUS (C=I): robusto ao acoplamento ---
    phis_raw_deg = zeros(1, L);
    A = zeros(M, L);
    for l = 1:L
        [~, phi_kw, ~] = doa_kw_uca(X_cell{l}, q.', radius, lambda, beta_uca);
        phis_raw_deg(l) = phi_kw(1);
        A(:,l) = utils.steering_vec_uca(M, radius, lambda, 90, phis_raw_deg(l));
    end

    % --- 2. Pool: UMA C das L assinaturas + direcoes (KW cru) ---
    [C_hat, ~] = estimate_C_circulant_pool(B_orig, A, M, selfcal_lambda, selfcal_prior);

    phis_hat_deg = phis_raw_deg;     % direcoes correntes (KW cru no inicio)

    % residuo observavel (SEM oracle) do estado corrente
    res_cur = pool_residual_local(B_orig, A, C_hat);
    history.res_per_iter = res_cur;          % log: residuo por iteracao
    % historico de direcoes por iteracao de refino (L x (n_refine+1)):
    % coluna 1 = KW cru (it 0); colunas seguintes = estado apos cada passo.
    history.phis_per_iter = phis_hat_deg(:);

    % --- 3. Refino conjunto C <-> direcoes, estilo Zhang ---
    %   A cada passo: recompensa -> KW refina direcoes -> re-estima C.
    %   use_safeguard=true : so' ACEITA o passo se o residuo observavel
    %     ||b_l - beta_l*C*a(phi_l)|| diminuir; senao REVERTE e para.
    %   use_safeguard=false: aplica todos os n_refine passos cegamente
    %     (para estudar a importancia do criterio de parada).
    history.stopped_at_iter = NaN;
    for it = 1:n_refine
        D = safe_inv_local(C_hat);
        Aref   = zeros(M, L);
        phiref = zeros(1, L);
        for l = 1:L
            Y = D * X_cell{l};
            [~, phi_kw, ~] = doa_kw_uca(Y, q.', radius, lambda, beta_uca);
            phiref(l) = phi_kw(1);
            Aref(:,l) = utils.steering_vec_uca(M, radius, lambda, 90, phiref(l));
        end
        [C_try, ~] = estimate_C_circulant_pool(B_orig, Aref, M, selfcal_lambda, selfcal_prior);
        res_try = pool_residual_local(B_orig, Aref, C_try);

        if ~use_safeguard
            % SEM salvaguarda: aceita sempre (refino cego)
            C_hat        = C_try;
            phis_hat_deg = phiref;
            res_cur      = res_try;
            history.res_per_iter(end+1) = res_cur;
            history.phis_per_iter(:,end+1) = phis_hat_deg(:);
        elseif res_try < res_cur * (1 - 1e-4)    % melhora real -> ACEITA
            C_hat        = C_try;
            phis_hat_deg = phiref;
            res_cur      = res_try;
            history.res_per_iter(end+1) = res_cur;
            history.phis_per_iter(:,end+1) = phis_hat_deg(:);
        else                                     % nao melhora -> REVERTE e para
            history.res_per_iter(end+1) = res_cur;   % registra o que ficou
            history.phis_per_iter(:,end+1) = phis_hat_deg(:);  % estado mantido
            history.stopped_at_iter = it;
            break;
        end
    end

    % --- log (sem oracle) ---
    history.phis_raw = phis_raw_deg;
    history.c_col    = C_hat(:,1);
    if ~isempty(C_true_for_log)
        history.errF = norm(C_hat - C_true_for_log,'fro') / max(norm(C_true_for_log,'fro'), eps);
    else
        history.errF = NaN;
    end
end

function D = safe_inv_local(C)
    M = size(C,1);
    rc = rcond(C);
    if isfinite(rc) && rc > 1e-12
        D = inv(C); if all(isfinite(D(:))), return; end
    end
    mu = 1e-6*(norm(C,'fro')^2/M + eps);
    D = (C'*C + mu*eye(M)) \ C';
    if ~all(isfinite(D(:))), D = eye(M); end
end

function r = pool_residual_local(B, A, C)
% Residuo observavel agregado (SEM oracle): quao bem C*a(phi_l) explica b_l,
% somado sobre as L transmissoes, com escala beta_l de melhor ajuste.
%   r = sum_l ||b_l - beta_l*C*a_l||^2 / sum_l ||b_l||^2,
%   beta_l = (C a_l)^H b_l / ||C a_l||^2
    L = size(B,2);
    num = 0; den = 0;
    for l = 1:L
        Ca = C * A(:,l);
        nb2 = real(B(:,l)'*B(:,l));
        nCa2 = real(Ca'*Ca) + eps;
        beta = (Ca' * B(:,l)) / nCa2;
        e = B(:,l) - beta*Ca;
        num = num + real(e'*e);
        den = den + nb2;
    end
    r = num / (den + eps);
end
