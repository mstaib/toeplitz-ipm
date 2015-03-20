function [ value, grad, hess ] = barrier(v)
    R = levinson_durbin(v);
    if nargout >=0
        value = barrier_value(v,R);
        %value = -sum(log(v));
    end
    if nargout >= 2
        grad = barrier_grad(v, R);
        %grad = -1./v;
    end
    if nargout >= 3
        hess = barrier_hess(v, R);
        %hess = diag(1./(v.^2));
    end
end

