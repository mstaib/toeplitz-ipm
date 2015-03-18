function [u, obj] = non_symm_primal_dual_ipm(c, A, b, x_feas, beta, gamma) 

v = phase1_central_path(A, b, x_feas, beta, gamma);
t = 1; %how to get this???

nu = length(v);
u = v;

obj = [c' * u];

for k=1:1000
    R = levinson_durbin(u);
    grad_u = barrier_grad(u, R);
    hess_u = barrier_hess(u, R);

    %% stopping criterion
    %[d, y] = delta(c, A, b, grad_u, hess_u, 0);
    %lambda = sqrt(d' * hess_u * d);
    %if lambda <= 2 * beta
    %    break
    %end

    % z_k = p_tau(u)
    [x, s, y] = p_tau(u, t, -1, grad_u, hess_u, c, A, b, gamma);

    % t = t(z_k)
    t = nu / (s' * x);

    % u = sigma_tau(z_k)
    u = sigma_tau(x, t, c, A, b, beta);

    obj = [obj c'*u];
end
