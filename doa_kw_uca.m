function [theta_deg, phi_deg, Bhat] = doa_kw_uca(X, Y, r, lambda, beta, opts)
% DOA_KW_UCA  Estimação 2D de DoA em UCA com formas de onda conhecidas
%
% Entradas:
%   X      : MxN, amostras do array (cada linha = sensor do UCA)
%   Y      : KxN, formas de onda conhecidas (cada linha = y_k(t) sincronizado)
%   r      : raio do UCA (em metros ou na mesma unidade de lambda)
%   lambda : comprimento de onda
%   beta   : 1xM, ângulos de cada sensor no círculo (rad). Use beta(m)=2*pi*(m-1)/M.
%   opts   : struct opcional:
%            .useUnwrap (default=false) – tenta desembrulhar fases (se r/lambda > 1/4)
%            .deltas    (default=1:floor(M/2)) – lista de saltos Δ usados
%            .regRyy    (default=1e-6) – regularização de Ryy (Tikhonov)
%
% Saídas:
%   theta_deg : 1xK elevação (graus, 0° = eixo z, 90° = plano do UCA)
%   phi_deg   : 1xK azimute (graus, -180..180)
%   Bhat      : MxK matriz de assinatura espacial estimada (B̂ = Rxy Ryy^{-1})
%


    arguments
        X {mustBeNumeric, mustBeNonempty}
        Y {mustBeNumeric, mustBeNonempty}
        r (1,1) double {mustBePositive}
        lambda (1,1) double {mustBePositive}
        beta (:,1) double
        opts.useUnwrap (1,1) logical = false
        opts.deltas = []
        opts.regRyy (1,1) double = 1e-6
    end

    [M, N] = size(X);
    [K, Ny] = size(Y);
    if Ny ~= N
        error('Y deve ter o mesmo número de snapshots que X (N).');
    end
    if numel(beta) ~= M
        error('beta deve ter dimensão M (um ângulo por sensor).');
    end

    % ----- 1) Estima B̂ = Rxy * (Ryy + μI)^{-1}  (Eq. 8)
    Rxy = (X * Y') / N;
    Ryy = (Y * Y') / N;
    if opts.regRyy > 0
        Ryy = Ryy + opts.regRyy * eye(K);
    end
    Bhat = Rxy / Ryy;  % MxK
    Bhat = Bhat ./ abs(Bhat(1,:));   % força referência de fase na antena 1

    % ----- 2) Preparar baselines Δ e matriz C empilhada (Eqs. 12-18)
    if isempty(opts.deltas)
        S = floor(M/2);
        deltas = 1:S;
    else
        deltas = opts.deltas(:).';
    end

    C_all = [];
    for Delta = deltas
        C_D = build_C_block(beta, Delta);  % M x 2 (com rotação circular)
        C_all = [C_all; C_D];              %#ok<AGROW>
    end
    % C_all tem dimensão (M*|deltas|) x 2

    % ----- 3) Estimar DoA por fonte k usando as fases de b̂_k  (Eqs. 11, 19, 20)
    theta_deg = zeros(1,K);
    phi_deg   = zeros(1,K);

    for k = 1:K
        bk = Bhat(:,k);  % assinatura espacial ~ a(θk, φk)

        eta_stack = [];
        for Delta = deltas
            % vetor eta_k(Δ) com M entradas (Eq. 12)
            eta_D = zeros(M,1);
            for m = 1:M
                p = m;
                q = wrap_index(m + Delta, M);  % índice cíclico
                ph = angle( bk(p) * conj(bk(q)) );
                eta_D(m) = ph;
            end
            if opts.useUnwrap
                % desembrulha em torno do círculo
                eta_D = unwrap_circular_phases(eta_D);
            end
            eta_stack = [eta_stack; eta_D]; %#ok<AGROW>
        end

        % LS: ĝ = (C^T C)^{-1} C^T eta   (Eq. 19)
        ghat = (C_all' * C_all) \ (C_all' * eta_stack);

        % Recupera (θ, φ) (Eq. 20)
        norm_g = hypot(ghat(1), ghat(2));
        sin_theta = (lambda/(2*pi*r)) * norm_g;
        % Clampeia por segurança numérica
        sin_theta = max(-1, min(1, sin_theta));
        theta = asin(sin_theta);

        phi = atan2(ghat(2), ghat(1));

        theta_deg(k) = rad2deg(theta);
        phi_deg(k)   = wrapTo180(rad2deg(phi));
    end
end

% ======================= Helpers =========================

function C = build_C_block(beta, Delta)
% Constrói o bloco C(Δ) de dimensão Mx2 (Eq. 13),
% considerando os pares (m, m+Δ) com índice cíclico.
    M = numel(beta);
    C = zeros(M,2);
    for m = 1:M
        p = m;
        q = wrap_index(m + Delta, M); % circular
        C(m,1) = cos(beta(p)) - cos(beta(q));
        C(m,2) = sin(beta(p)) - sin(beta(q));
    end
end

function idx = wrap_index(i, M)
% Índice circular em {1,...,M}
    idx = mod(i-1, M) + 1;
end

function y = unwrap_circular_phases(x)
% Desembrulha um vetor de M fases dispostas ao longo do círculo do UCA.
% Faz unwrap na "ordem física" dos sensores e depois reancora para [-pi,pi].
    y = unwrap(x);
    % normaliza para [-pi, pi]
    y = mod(y + pi, 2*pi) - pi;
end
