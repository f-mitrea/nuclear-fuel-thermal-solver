function[F] = compute_pellet_residual(T,radius,dr,q_vol,T_wall,k_UO2,dk_dT)

N = length(T);
F = zeros(N,1);

F(1) = 4*k_UO2(T(1)) * (T(2) - T(1))/(dr^2) + q_vol;

for i = 2:N-1
    dT_dr = (T(i+1)-T(i-1))/(2*dr);
    d2T_dr2 = (T(i+1) - 2*T(i) + T(i-1))/(dr^2);
    F(i) = k_UO2(T(i))*d2T_dr2 + k_UO2(T(i))/radius(i)*dT_dr + dk_dT(T(i))*(dT_dr^2) + q_vol;
end

F(N) = T(N) - T_wall;

end