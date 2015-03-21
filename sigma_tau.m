function u = sigma_tau(x, tau, c, A, b, beta)

u = x;

last_two_lambdas = [];
iters = 1;

while 1
    R = levinson_durbin(u);
    grad_u = barrier_grad(u, R);
    hess_u = barrier_hess(u, R);
    [d_u, y_u] = delta(c, A, b, grad_u, hess_u, tau);
    %min(eig(toeplitz2(u)));

    lambda = sqrt(d_u' * hess_u * d_u)
    if lambda <= beta
        break
    end

    if length(last_two_lambdas) >= 2
        if sum(abs(last_two_lambdas - lambda) < 1e-3) || lambda > last_two_lambdas(2)
            break
        end
        last_two_lambdas = [last_two_lambdas(2) lambda];
    else
        last_two_lambdas = [last_two_lambdas lambda];
    end


    u = u + d_u / (1 + lambda);
    iters = iters + 1;
end
