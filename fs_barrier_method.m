function [ x ] = fs_barrier_method(c, A, b, x_0, epsilon, mu, barrier, theta)
%BARRIER_METHOD Summary of this function goes here
    n = size(x_0,1);
    t = 1;
    x = x_0;
    if ~exist('eps','var') || isempty(mu)
        epsilon = 1e-3;
    end
    if ~exist('mu','var') || isempty(mu)
        mu=20;
    end
    function value = f_val(c,x)
        value = c'*x + barrier(x);
    end

    function [grad, hess] = f_derivs(c,x)
        [~,grad,hess] = barrier(x);
        grad = grad + c;
    end
    history = [];
    iters = 0;
    while true
        [x_star, nu_star, steps] = fs_newton(A,b,t*c,x, @f_val, @f_derivs);
        x = x_star;
        history = [history [steps; n/t]];
        if theta/t < epsilon
            break
        end
        iters = iters + 1;
        t = mu*t;
        fprintf('outer iteration\n')
    end
    fprintf('Converged.\n  outer iterations = %d, newton steps = %d\n  optval = %f, duality gap = %f\n', iters, sum(history(1,:)), c'*x_star, n/t)
    nu_star = nu_star/t;
end

