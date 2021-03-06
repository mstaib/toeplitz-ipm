function [x_star, nu_star, nt_steps] = fs_newton(A, b, c, x_0, f_val, f_derivs, alpha, beta)
%FS_NEWTON Feasible start Newton's method
    x = x_0;
    nt_steps = 0;
    if ~exist('beta','var') || isempty(beta)
        beta=0.5;
    end
    if ~exist('alpha','var') || isempty(alpha)
        alpha=0.01;
    end
    epsilon = 1e-6;
    while true
       [grad, hess] = f_derivs(c,x);
       %fprintf('hessian condition number: %f\n',cond(hess))
       %fprintf('number of small hessian eigenvals: %d\n',sum(eig(hess)<=1e-3))
       % compute newton step
       [nt_dir, w] = solve_kkt(hess,grad,A);
       lambda_sq = -nt_dir'*grad;
       lambda_sq
       if lambda_sq/2 <= epsilon
           break;
       end
       % line search
       t = 1;
       while det(toeplitz2(x)) < 0
           t = beta*t;
       end
       val = f_val(c,x);
       while (f_val(c, x + t*nt_dir) > val+ alpha*t*grad'*nt_dir)
           t = beta*t;
       end
       % update
       x = x + t*nt_dir;
       nt_steps = nt_steps + 1;
       fprintf('inner iteration\n')
    end
    x_star = x;
    nu_star = w;
    % check KKT
    %assert(norm(A*x_star -b) <= 1e-10);
    %assert(norm(f_grad(c,x_star) + A'*nu_star) <= 1e-2);
end

function [nt_dir, w] = solve_kkt(H,g,A)
    [m,n] = size(A);
    H = H + A'*A;
    kkt_A = [H A';A zeros(m,m)];
    kkt_b = [-g;zeros(m,1)];
    x = kkt_A\kkt_b;
    nt_dir = x(1:length(g));
    w = x(length(g)+1:end);
end

function [x_1, x_2] = block_elim(hess, A, g)
    [m,n] = size(A);
    hess = hess + A'*eye(m)*100*A;
    cond(hess)
    S = -A*(hess\A');
    b = A*(hess\g);
    x_2 = S\b;
    x_1 = hess\(-A'*x_2 - g);
end