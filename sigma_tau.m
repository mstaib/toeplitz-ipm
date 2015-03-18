function u = sigma_tau(x, tau, c, A, b, beta)

u = x;
while 1
    R = levinson_durbin(u);
    grad_u = barrier_grad(u, R);
    hess_u = barrier_hess(u, R);
    [d_u, y_u] = delta(c, A, b, grad_u, hess_u, tau);

    lambda = sqrt(d_u' * hess_u * d_u)
    if lambda <= beta
        break
    end

    u = u + d_u / (1 + lambda);
end
