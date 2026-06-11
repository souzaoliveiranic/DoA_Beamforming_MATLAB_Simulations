function d_deg = spherical_angular_distance(theta1_deg, phi1_deg, theta2_deg, phi2_deg)
    % SPHERICAL_ANGULAR_DISTANCE  Distancia angular esferica (geodesica)
    % entre dois pontos definidos por (theta, phi) em graus.
    %
    %   theta  = angulo de elevacao a partir do eixo +z (zenith=0, plano XY=90)
    %   phi    = angulo de azimute no plano XY
    %
    %   Vetor unitario: u = [sin(theta)cos(phi), sin(theta)sin(phi), cos(theta)]
    %   Distancia: d = acos(u1' * u2)
    %
    %   Resultado em graus, no intervalo [0, 180].

    t1 = deg2rad(theta1_deg);  p1 = deg2rad(phi1_deg);
    t2 = deg2rad(theta2_deg);  p2 = deg2rad(phi2_deg);

    u1 = [sin(t1)*cos(p1); sin(t1)*sin(p1); cos(t1)];
    u2 = [sin(t2)*cos(p2); sin(t2)*sin(p2); cos(t2)];

    c = max(min(u1.' * u2, 1), -1);   % clamp para [-1,1] (estabilidade numerica)
    d_deg = rad2deg(acos(c));
end
