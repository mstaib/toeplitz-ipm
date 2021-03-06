clear all; close all;
% set up
m = 10;
n = 20;
cvx_status = 'Infeasible';
while strcmp(cvx_status,'Infeasible')
    seed=randi(10000);
    rng(seed);
    A = rand(m,n);
    x = rand(n,1);
    b = A*x;
    c = rand(n,1);
    fprintf('Creating a feasible starting point...\n')
    cvx_begin quiet
        variable X(n,n) Toeplitz Semidefinite
        minimize 0
        subject to
            A*(X(1,:)'.*[1/2; ones(n-1,1)]) == b
    cvx_end
end


x_0 = X(1,:)';
x_0(1) = x_0(1)/2;
fprintf('Solving with the barrier method.\n')
mu=20;
tic
x_bm = fs_barrier_method(c, A, b, x_0, 1e-4, mu, @barrier, n+1);
toc

fprintf('Solving with cvx.\n')
tic
cvx_begin quiet
    variable X(n,n) Toeplitz Semidefinite
    minimize (c'*X(1,:)')
    subject to
        A*X(1,:)' == b
cvx_end
toc
x_cvx = X(1,:)';

norm(x_cvx - x_bm)