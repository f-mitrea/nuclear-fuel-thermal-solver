clc,clear,close all

addpath('C:\Users\mfale\Desktop\nuclear-fuel-thermal-solver\matlab_prototype\common');

syms T
k_UO2 = input('Enter the UO2 thermal conductivity [W/(m*K)]: ');
if isempty(k_UO2)
    k_UO2 = 1/(0.0375 + 2.165e-4*T);
end
k_sym = sym(k_UO2);
dk_sym = diff(k_sym,T);
d2k_sym = diff(diff(k_sym,T));
k_UO2 = matlabFunction(k_sym,'Vars',T);
dk_dT = matlabFunction(dk_sym,'Vars',T);
d2k_dT2 = matlabFunction(d2k_sym,'Vars',T);

T_wall = input('Enter the outer boundary temperature [K]: ');
if isempty(T_wall)
    T_wall = 600;
end

q_vol = input('Enter the volumetric heat generation [W/m^3]: ');
if isempty(q_vol)
    q_vol = 150e6;
end

R_pellet = input('Enter the fuel pellet radius [m]: ');
if isempty(R_pellet)
    R_pellet = 4e-3;
end

n_nodes = input('Enter the number of radial nodes: ');
if isempty(n_nodes)
    n_nodes = 100;
end

tol = input('Enter the Newton convergence tolerance: ');
if isempty(tol)
    tol = 1e-8;
end

max_iter = input('Enter the Newton number of iterations: ');
if isempty(max_iter)
    max_iter = 7;
end

dr = R_pellet/(n_nodes-1);
radius = zeros(n_nodes,1);

for i = 1:n_nodes
    radius(i) = (i-1)*dr;
end

T = zeros(n_nodes,1);

for i = 1:n_nodes
    T(i) = T_wall + q_vol*(R_pellet^2 - radius(i)^2)/(4*k_UO2(T_wall));
end

[T] = solve_pellet_newton(T,radius,dr,tol,max_iter,k_UO2,dk_dT,d2k_dT2,q_vol,T_wall);

figure
plot(radius*1e3, T, 'o')
xlabel('r [mm]')
ylabel('T [K]')
grid on



