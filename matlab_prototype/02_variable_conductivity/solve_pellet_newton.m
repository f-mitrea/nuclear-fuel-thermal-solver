function[T_new] = solve_pellet_newton(T_old,radius,dr,tol,max_iter,k_UO2,dk_dT,d2k_dT2,q_vol,T_wall)

n_iter = 0;
err = 1;

while norm(err,inf) > tol && n_iter < max_iter
    [J] = assemble_pellet_jacobian(T_old,radius,dr,k_UO2,dk_dT,d2k_dT2);
    [F] = compute_pellet_residual(T_old,radius,dr,q_vol,T_wall,k_UO2,dk_dT);
    b = -F;
    [delta_T] = solve_tridiagonal_system(J,b);
    T_new = T_old + delta_T;
    err = T_new - T_old;
    T_old = T_new;
    n_iter = n_iter + 1;
end

end
