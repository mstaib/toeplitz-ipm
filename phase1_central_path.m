function v = phase1_central_path(A, b, x_feas, beta, gamma) 

nu = length(x_feas);

v = x_feas;
tau = 1;

R = levinson_durbin(v);
grad_v = barrier_grad(v, R);
c = -grad_v;

while 1
    R = levinson_durbin(v);
    grad_v = barrier_grad(v, R);
    hess_v = barrier_hess(v, R);

    % stopping criterion
    [d, y] = delta(c, A, b, grad_v, hess_v, 0);
    lambda = sqrt(d' * hess_v * d)
    if lambda <= 2 * beta
        break
    end

    % z_k = p_tau^+(v)
    [x, s, y] = p_tau(v, tau, 1, grad_v, hess_v, c, A, b, beta, gamma);

    % tau = t(z_k)
    tau = nu / (s' * x)

    % v = sigma_tau(z_k)
    % should be tau instead of 4/tau??
    v = sigma_tau(x, nu/tau, c, A, b, beta);
end
