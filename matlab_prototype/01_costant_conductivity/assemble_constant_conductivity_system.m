function[A,b] = assemble_constant_conductivity_system(p_coeff,q_coeff,r_coeff,dr,n_nodes,T_wall,q_vol,k_UO2)

A = zeros(n_nodes,n_nodes);
b = zeros(n_nodes,1);
gamma = q_vol*dr^2/k_UO2;

for i = 1:n_nodes
    if i > 1 && i < n_nodes
        b(i) = dr^2*r_coeff(i);
        A(i,i) = 2+dr^2*q_coeff(i);
        A(i,i+1) = -1+dr/2*p_coeff(i);
        A(i,i-1) = -1-dr/2*p_coeff(i);
    elseif i == 1
        b(i) = gamma;
        A(i,i) = 4;
        A(i,i+1) = -4;
    elseif i == n_nodes
        b(i) = dr^2*r_coeff(i)+(1-dr/2*p_coeff(i))*T_wall;
        A(i,i) = 2+dr^2*q_coeff(i);
        A(i,i-1) = -1-dr/2*p_coeff(i);
    end
end

end