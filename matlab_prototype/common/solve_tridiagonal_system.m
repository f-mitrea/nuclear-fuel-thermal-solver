function[Y] = solve_tridiagonal_system(A,B)

N = size(A,1);
X = zeros(N,1);
Y = zeros(N,1);
u = zeros(N,1);
v = zeros(N,1);
alphas = zeros(N,1);

for i = 1:N-1
    v(i) = A(i,i+1);
end

u(1) = A(1,1);

for j = 2:N
    alphas(j) = A(j,j-1)/u(j-1);
    u(j) = A(j,j) - alphas(j)*v(j-1);
end

%L = eye(N) + diag(alphas(2:N),-1);
%U = diag(u) + diag(v(1:N-1),1);

X(1) = B(1);

for k = 2:N
    X(k) = B(k) - alphas(k)*X(k-1);
end

Y(N) = X(N)/u(N);

for l = N-1:-1:1
    Y(l) = (X(l)-v(l)*Y(l+1))/u(l);
end

end