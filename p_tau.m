function [x, s, y] = p_tau(u, tau, multiplier, grad_u, hess_u, c, A, b, beta, gamma)

[d, y] = delta(c, A, b, grad_u, hess_u, tau);

xk = u - d;
sk = c - 1/tau * A' * y;
yk = 1/tau * y;

w = sqrt(tau) * u;
R = levinson_durbin(w);
hess_w = barrier_hess(w, R);

[del_x, del_s, del_y] = delta_z(c, A, b, hess_w, sk);

alpha = get_stepsize(xk, sk, yk, del_x, del_s, del_y, beta, gamma);
x = xk + multiplier * alpha * del_x;
s = sk + multiplier * alpha * del_s;
y = yk + multiplier * alpha * del_y;

end
