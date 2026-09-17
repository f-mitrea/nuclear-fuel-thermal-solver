clc,clear,close all

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));

k_UO2 = input('Enter the UO2 thermal conductivity [W/(m*K)]: ');
if isempty(k_UO2)
    k_UO2 = 2.3;
end

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


dr = R_pellet/n_nodes;

radius = zeros(n_nodes,1);
for i = 1:n_nodes
    radius(i) = (i-1)*dr;
end


p_coeff = -1./radius;
p_coeff(1) = 0;

r_coeff = q_vol/k_UO2 * ones(n_nodes,1);

q_coeff = zeros(n_nodes,1);

[A,b] = assemble_constant_conductivity_system(p_coeff,q_coeff,r_coeff,dr,n_nodes,T_wall,q_vol,k_UO2);
[T_num] = solve_tridiagonal_system(A,b);

radius(n_nodes+1) = R_pellet;
T_num(n_nodes+1) = T_wall;

T_exact = T_wall + q_vol*(R_pellet^2-radius.^2)./(4*k_UO2);

err = abs(T_num - T_exact);

figure
plot(radius*1e3, T_num, 'o', radius*1e3, T_exact, '-')
xlabel('Radius [mm]');
ylabel('Temperature [K]')
legend('Numerical','Analytical')
grid on

figure
plot(radius*1e3, err, '-')
xlabel('Radius [mm]');
ylabel('T_{num} - T_{exact} [K]')
grid on





