function [C_hat, c_hat, alpha_hat, phi_hat_deg, history] = ...
    estimate_C_selfcal(X, q, M, radius, lambda, doa_method, ...
                       phi_grid_deg, max_iter, tol_phi, tol_C, C_true_for_log, ...
                       damping, init_mode)
% ESTIMATE_C_SELFCAL  Self-calibration alternante: estima conjuntamente
%   DoA e matriz de acoplamento. Aceita 4 metodos de DoA: 'DAS', 'CAPON',
%   'MUSIC', 'KW'.
%
% Entradas:
%   X            : MxN, dados recebidos
%   q            : Nx1 ou 1xN, forma de onda do sinal util (so usado p/ KW)
%   M            : numero de elementos
%   radius       : raio do UCA (em metros)
%   lambda       : comprimento de onda (em metros)
%   doa_method   : 'DAS', 'CAPON', 'MUSIC' ou 'KW'  (case-insensitive)
%   phi_grid_deg : grid de busca para o DoA  (padrao -180:0.5:180)
%   max_iter     : max. iteracoes alternantes  (padrao 15)
%   tol_phi      : tolerancia em |phi_t - phi_{t-1}| (graus)  (padrao 0.01)
%   tol_C        : tolerancia em ||C - C_prev||_F / ||C_prev||_F  (padrao 1e-8)
%   C_true_for_log : (opcional) C_true para logging do erro Frobenius
%                    relativo a cada iteracao. Se vazio, history.err_F = NaN.
%   damping      : (opcional, default 1.0) fator de subrelaxacao em [0,1].
%                  C_hat <- (1-damping)*C_prev + damping*C_new
%                  damping=1 => sem damping (default classico).
%                  damping=0.5 => mistura 50-50, suaviza oscilacoes.
%   init_mode    : (opcional, default 'identity'). Pode ser:
%                  - 'identity' : C_0 = I_M (init classico)
%                  - 'kw'       : estima phi via doa_kw_uca SEM acoplamento
%                                 e gera C_0 via estimate_C_circulant_uca
%                                 (recomendado quando ha interferentes fortes)
%                  - matriz MxM : usa diretamente como C_0 (ex.: matriz da
%                                 calibracao offline)
%
% Saidas:
%   C_hat       : MxM, matriz de acoplamento estimada
%   c_hat       : K x 1, coeficientes circulantes unicos (K = floor(M/2))
%   alpha_hat   : escalar complexo (ganho)
%   phi_hat_deg : DoA estimada na ultima iteracao
%   history     : struct com tracking por iteracao:
%                 .phi_per_iter   - vetor de DoAs estimadas
%                 .err_F_per_iter - vetor de erros Frobenius (se C_true dado)
%                 .delta_C        - vetor de mudancas em ||C - C_prev||_F
%                 .delta_phi      - vetor de |phi_t - phi_{t-1}|
%                 .n_iter         - num. iteracoes ate convergencia
%
% Algoritmo (alternante):
%   Init: C = I_M, phi inicial = DoA(X)
%   Repeat:
%     1) D = inv(C);  Y = D*X
%     2) phi <- DoA(Y) usando metodo selecionado:
%          - 'DAS', 'CAPON', 'MUSIC': grid search no espectro espacial de Y
%          - 'KW': usa doa_kw_uca (forma fechada 2D, baseada em fases
%                  diferenciais) sobre Y. Mais rapido e mais preciso que
%                  grid search; e' o mesmo estimador que o usuario usa no
%                  loop principal, garantindo comparacao justa entre
%                  Self-cal-KW e KW LS 1-dir offline.
%     3) C <- estimate_C_circulant_uca(X q* / q'q, a(phi), M)
%        OBS: o LS de C usa SEMPRE b_hat ORIGINAL (X*q*/q'q), nao Y*q*/q'q.
%        Isto porque o modelo eh b = alpha * C * a + ruido; aplicar D antes
%        do KW so faz sentido para a estimacao de DoA, nao para a de C.
%   Until convergencia em phi e C.

    if nargin < 6 || isempty(doa_method),    doa_method = 'KW';            end
    if nargin < 7 || isempty(phi_grid_deg),  phi_grid_deg = -180:0.5:180;  end
    if nargin < 8 || isempty(max_iter),      max_iter = 15;                end
    if nargin < 9 || isempty(tol_phi),       tol_phi  = 0.01;              end
    if nargin < 10 || isempty(tol_C),        tol_C    = 1e-8;              end
    if nargin < 11,                          C_true_for_log = [];          end
    if nargin < 12 || isempty(damping),      damping = 1.0;                end
    if nargin < 13 || isempty(init_mode),    init_mode = 'identity';       end

    doa_method = upper(doa_method);
    q = q(:);
    qHq = q' * q;

    % b_hat sem compensacao: usado tanto para DoA-KW quanto para LS de C
    b_hat_orig = X * conj(q) / qHq;

    % Pre-calcula o dicionario de steering vectors (acelera muito o grid search)
    nPhi = numel(phi_grid_deg);
    A_dict = zeros(M, nPhi);
    for ip = 1:nPhi
        A_dict(:, ip) = utils.steering_vec_uca(M, radius, lambda, 90, phi_grid_deg(ip));
    end
    norm_a_sq = real(sum(conj(A_dict).*A_dict, 1)).';   % nPhi x 1

    % --- Inicializacao --------------------------------------------------
    if isnumeric(init_mode) && all(size(init_mode) == [M M])
        % init_mode e' uma matriz MxM pre-computada (ex.: cal. offline)
        C_hat = init_mode;
    else
        init_mode = lower(init_mode);
        switch init_mode
            case 'identity'
                C_hat = eye(M);
            case 'kw'
                % Estima phi inicial via doa_kw_uca direto sobre X (sem D),
                % depois gera C_0 via LS circulante. Como doa_kw_uca usa
                % apenas fases (relativamente robustas a acoplamento moderado),
                % isto da' um ponto de partida muito melhor que C=I quando
                % ha interferentes fortes.
                beta_uca_init = 2*pi*(0:M-1).'/M;
                opts_kw_init  = struct('regRyy', 1e-6);
                [~, phi0_deg, ~] = doa_kw_uca(X, q.', radius, lambda, ...
                                              beta_uca_init);
                a0 = utils.steering_vec_uca(M, radius, lambda, 90, phi0_deg(1));
                [C_hat, ~, ~, ~] = estimate_C_circulant_uca(b_hat_orig, a0, M);
            otherwise
                error('estimate_C_selfcal:init', ...
                      ['init_mode invalido: %s. Use ''identity'', ''kw'' ' ...
                       'ou uma matriz MxM.'], init_mode);
        end
    end
    phi_prev = NaN;
    C_prev   = C_hat;

    % --- Tracking de melhor C visto (para "early stopping" se divergir) -
    if ~isempty(C_true_for_log)
        best_err_F = norm(C_hat - C_true_for_log, 'fro') / ...
                     max(norm(C_true_for_log,'fro'), eps);
        best_C     = C_hat;
    else
        best_err_F = inf;
        best_C     = C_hat;
    end

    history.phi_per_iter   = zeros(max_iter, 1);
    history.err_F_per_iter = nan(max_iter, 1);
    history.delta_C        = nan(max_iter, 1);
    history.delta_phi      = nan(max_iter, 1);
    history.n_iter         = 0;

    if ~isempty(C_true_for_log)
        norm_Ctrue = norm(C_true_for_log, 'fro');
    end

    n_iter = 0;
    for it = 1:max_iter
        % --- 1. Aplica D = inv(C) e prepara dados compensados p/ DoA ----
        D = inv(C_hat);
        Y = D * X;                            % dados compensados (todos os metodos usam)

        % --- 2. Estima DoA usando o metodo selecionado ------------------
        switch doa_method
            case 'KW'
                % Usa o estimador KW 2D em forma fechada do usuario (doa_kw_uca).
                % Diferente dos outros metodos, aqui obtemos phi continuo
                % (sem viés de quantizacao do grid).
                Y_kw = q.';                       % 1 x N (1 forma de onda)
                beta_uca = 2*pi*(0:M-1).'/M;
                opts_kw = struct('regRyy', 1e-6); % mesmo default da utils
                [~, phi_kw_deg, ~] = doa_kw_uca(Y, Y_kw, radius, lambda, ...
                                                beta_uca);
                phi_hat_deg_kw = phi_kw_deg(1);
                % Marca para o codigo abaixo nao usar A_dict(ip_best)
                ip_best = NaN;

            case 'DAS'
                R_y = (Y * Y') / size(Y,2);
                scores = real(sum(conj(A_dict) .* (R_y * A_dict), 1)).' ./ ...
                         max(norm_a_sq, eps);
                [~, ip_best] = max(scores);

            case 'CAPON'
                R_y = (Y * Y') / size(Y,2);
                R_y = R_y + 1e-6 * trace(R_y) / M * eye(M);  % regulariza
                Rinv = R_y \ eye(M);
                den = real(sum(conj(A_dict) .* (Rinv * A_dict), 1)).';
                den = max(den, eps);
                scores = 1 ./ den;
                [~, ip_best] = max(scores);

            case 'MUSIC'
                R_y = (Y * Y') / size(Y,2);
                [V, D_eig] = eig(R_y);
                [~, idx_sort] = sort(real(diag(D_eig)), 'descend');
                V = V(:, idx_sort);
                % Assume 1 fonte util + 1 interferente => P=2
                P_sources = 2;
                E_n = V(:, P_sources+1:end);
                proj = E_n' * A_dict;          % (M-P) x nPhi
                den = real(sum(conj(proj).*proj, 1)).';
                den = max(den, eps);
                scores = 1 ./ den;
                [~, ip_best] = max(scores);

            otherwise
                error('estimate_C_selfcal:method', ...
                      'Metodo desconhecido: %s. Use DAS|CAPON|MUSIC|KW.', doa_method);
        end

        % --- 2.1 Recupera phi_hat e a_hat -------------------------------
        if strcmp(doa_method, 'KW')
            phi_hat_deg = phi_hat_deg_kw;
            a_hat = utils.steering_vec_uca(M, radius, lambda, 90, phi_hat_deg);
        else
            phi_hat_deg = phi_grid_deg(ip_best);
            a_hat = A_dict(:, ip_best);
        end

        % --- 3. Re-estima C com b_hat_orig e a(phi_hat) -----------------
        [C_new, c_hat, alpha_hat, ~] = estimate_C_circulant_uca(b_hat_orig, a_hat, M);

        % --- 3.1 Damping (subrelaxacao): suaviza oscilacoes -------------
        % C_t = (1-gamma)*C_{t-1} + gamma*C_new
        % Damping=1 => sem damping (default classico).
        % Damping<1 => mistura com iteracao anterior, reduz oscilacao.
        if damping ~= 1.0
            C_new = (1 - damping) * C_prev + damping * C_new;
            % Recalcula c_hat consistente com C_new amortecido
            c_hat = [C_new(1, 2:M/2).' ; C_new(1, M/2+1)];
        end

        % --- Logging por iteracao ---------------------------------------
        history.phi_per_iter(it) = phi_hat_deg;
        cur_err_F = nan;
        if ~isempty(C_true_for_log)
            cur_err_F = norm(C_new - C_true_for_log, 'fro') / norm_Ctrue;
            history.err_F_per_iter(it) = cur_err_F;
            % Atualiza melhor C visto (criterio anti-divergencia)
            if cur_err_F < best_err_F
                best_err_F = cur_err_F;
                best_C     = C_new;
            end
        end
        d_C   = norm(C_new - C_prev, 'fro') / max(norm(C_prev,'fro'), eps);
        d_phi = abs(phi_hat_deg - phi_prev);
        history.delta_C(it)   = d_C;
        history.delta_phi(it) = d_phi;

        % --- Atualiza estado --------------------------------------------
        C_hat    = C_new;
        phi_prev = phi_hat_deg;
        C_prev   = C_new;
        n_iter   = it;

        % --- Criterio de parada (apos pelo menos 2 iter) ----------------
        if (it >= 2) && (d_phi < tol_phi) && (d_C < tol_C)
            break;
        end
    end

    history.n_iter = n_iter;
    % Trunca historico ate iteracao real
    history.phi_per_iter   = history.phi_per_iter(1:n_iter);
    history.err_F_per_iter = history.err_F_per_iter(1:n_iter);
    history.delta_C        = history.delta_C(1:n_iter);
    history.delta_phi      = history.delta_phi(1:n_iter);

    % --- Anti-divergencia: se o C final for pior que o melhor visto, ----
    %     retorna o melhor (so faz sentido se C_true_for_log foi dado).  --
    if ~isempty(C_true_for_log)
        final_err_F = norm(C_hat - C_true_for_log, 'fro') / norm_Ctrue;
        if final_err_F > 1.5 * best_err_F
            % Algoritmo divergiu apos atingir o minimo. Volta para best_C.
            C_hat = best_C;
            history.diverged = true;
            history.best_err_F = best_err_F;
        else
            history.diverged = false;
            history.best_err_F = final_err_F;
        end
    end
end
