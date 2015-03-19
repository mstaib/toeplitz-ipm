function [u, obj, uhist] = non_symm_primal_dual_ipm(c, A, b, x_feas, beta, gamma) 

v = x_feas;
v = phase1_central_path(A, b, x_feas, beta, gamma);
t = 1; %how to get this???

nu = length(v);
u = v;

%R = levinson_durbin(v);
%grad_v = barrier_grad(v, R);
%hess_v = barrier_hess(v, R);
%[d, y] = delta(c, A, b, grad_v, hess_v, 0.95);
%lambda = sqrt(d' * hess_v * d)

break


uhist = [u];
obj = [c' * u];

for k=1:10
    R = levinson_durbin(u);
    grad_u = barrier_grad(u, R);
    hess_u = barrier_hess(u, R);

    % z_k = p_tau(u)
    [x, s, y] = p_tau(u, t, -1, grad_u, hess_u, c, A, b, beta, gamma);

    % t = t(z_k)
    t = nu / (s' * x);

    % u = sigma_tau(z_k)
    u = sigma_tau(x, t, c, A, b, beta);

    uhist = [uhist u];
    obj = [obj c'*u];
end
