function a = steering_vec(M, d, lambda, theta_deg)
% Vetor de direção ULA (campo distante, banda estreita), ângulo a partir da BROADSIDE.
% a(m) = exp(-j*2π (m-1) d/λ * sin(theta)), m=1..M
    m = (0:M-1).';
    theta = deg2rad(theta_deg);
    a = exp(-1j*2*pi*(d/lambda)*m*sin(theta));
    a = a / norm(a);  % normalização numérica
end