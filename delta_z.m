function [del_x, del_s, del_y] = delta_z(c, A, b, hess_w, s);

[m, n] = size(A);

left_mat = [hess_w, eye(n), zeros(n,m); A, zeros(m,n), zeros(m,m); zeros(n,n), eye(n), A'];
right_vec = [s; zeros(m,1); zeros(n,1)];

del_z = left_mat\right_vec;
del_x = del_z(1:n);
del_s = del_z((n+1):(2*n));
del_y = del_z((2*n+1):end);
