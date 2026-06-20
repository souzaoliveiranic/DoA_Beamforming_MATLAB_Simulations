% Diagnostico rapido da estimacao de C: condicionamento vs raio
clc;
M = 8; fc = 500e6; c = 3e8; lambda = c/fc;
for r_lam = [0.1 0.2 0.25 0.5]
    r = r_lam*lambda;
    C_true = compute_Ctx_for_R(fc, M, r, 50);
    cu = C_true(1, 2:M/2+1);   % c1..c4 (unicos)
    fprintf('\n===== r = %.2f lambda  (k*r = %.3f rad = %.1f deg) =====\n', ...
        r_lam, 2*pi*r/lambda, rad2deg(2*pi*r/lambda));
    fprintf('  ||C_true||_F = %.3f   cond(C_true) = %.2e\n', norm(C_true,'fro'), cond(C_true));
    for k=1:M/2
        fprintf('  c%d = %+.4f %+.4fj  (|c%d|=%.4f)\n', k, real(cu(k)), imag(cu(k)), k, abs(cu(k)));
    end

    % --- Teste de recuperacao SEM ruido, direcoes VERDADEIRAS ---
    rng(1);
    P = 5; phis = -180 + 360*rand(1,P);
    A = zeros(M,P); B = zeros(M,P);
    for p=1:P
        a = utils.steering_vec_uca(M, r, lambda, 90, phis(p));
        A(:,p) = a; B(:,p) = C_true*a;     % alpha=1, sem ruido
    end
    [C_hat,~,~,~] = estimate_C_multidir(B, A, M, 50, 1e-12);
    errF_clean = norm(C_hat - C_true,'fro')/norm(C_true,'fro');
    fprintf('  [SEM ruido, dir verdadeiras] Frobenius rel = %.3e\n', errF_clean);

    % --- condicionamento do LS empilhado (mesma construcao do estimate_C_multidir) ---
    K = M/2; M_big = zeros(M*P,K);
    for p=1:P
        for i=1:M
            for k=1:K
                ip = mod(i-1+k,M)+1; im = mod(i-1-k,M)+1;
                if k==K, M_big((p-1)*M+i,k)=A(ip,p);
                else,    M_big((p-1)*M+i,k)=A(ip,p)+A(im,p); end
            end
        end
    end
    fprintf('  cond(M_big) do LS de c = %.3e\n', cond(M_big));
end

function Ctx = compute_Ctx_for_R(fc, M, R, Z0)
    c = 3e8; lambda = c/fc;
    mp = dipole; mp.Length = 0.5*lambda; mp.Width = 0.01*lambda;
    uca = circularArray; uca.Element = mp; uca.NumElements = M; uca.Radius = R;
    Sobj = sparameters(uca, fc);
    Z_matrix = s2z(Sobj.Parameters(:,:,1), Z0);
    denom = diag(Z_matrix) + Z0;
    Ctx = eye(M);
    for j = 1:M
        for i = 1:M
            if i ~= j, Ctx(i,j) = Z_matrix(i,j) / denom(j); end
        end
    end
end
