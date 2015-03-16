function [d, y] = delta(c, A, b, grad_u, hess_u, t)

[m, n] = size(A);

left_mat = [hess_u, -A'; A, zeros(m)];
right_vec = [-t*c - grad_u; zeros(m,1)];

res = left_mat\right_vec;
d = res(1:n);
y = res((n+1):end);
