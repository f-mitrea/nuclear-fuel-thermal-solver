function[J] = assemble_pellet_jacobian(T,radius,dr,k_UO2,dk_dT,d2k_dT2)

N = length(T);
J = zeros(N,N);

J(1,1) = (4*dk_dT(T(1))*(T(2)-T(1))-4*k_UO2(T(1)))/(dr^2);
J(1,2) = 4*k_UO2(T(1))/(dr^2);
J(N,N) = 1;

for i = 2:N-1
    dT_dr = (T(i+1)-T(i-1))/(2*dr);
    d2T_dr2 = (T(i+1) - 2*T(i) + T(i-1))/(dr^2);
    J(i,i-1) = k_UO2(T(i))/(dr^2) - k_UO2(T(i))/(2*dr*radius(i)) - dk_dT(T(i))*dT_dr/dr;
    J(i,i+1) = k_UO2(T(i))/(dr^2) + k_UO2(T(i))/(2*dr*radius(i)) + dk_dT(T(i))*dT_dr/dr;
    J(i,i) = -2*k_UO2(T(i))/(dr^2) + dk_dT(T(i))*d2T_dr2 + dk_dT(T(i))/radius(i)*dT_dr + d2k_dT2(T(i))*dT_dr^2;
end

end


