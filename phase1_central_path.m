function v = phase1_central_path(c, A, b, x_feas) 

v = x_feas;
tau = 1;

% z_k = p_tau^+(v)
R = levinson_durbin(v);
grad_v = barrier_grad(v, R);
hess_v = barrier_hess(v, R);

[d, y] = delta(c, A, b, grad_v, hess_v, tau);

xk = v - d;
sk = c - 1/tau * A' * y;
yk = 1/tau * y;

w = sqrt(tau) * v;
R = levinson_durbin(w);
hess_w = barrier_hess(w, R);

[del_x, del_s, del_y] = delta_z(c, A, b, hess_w, sk);

% tau = t(z_k)

% v = sigma_tau(z_k)
