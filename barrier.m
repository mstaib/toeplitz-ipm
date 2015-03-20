function [ value, grad, hess ] = barrier(v)
%BARRIER_DERIVS Summary of this function goes here
    R = levinson_durbin(v);
    if nargout >=1
        value = barrier_value(v,R);
        %value = -sum(log(v));
    end
    if nargout >= 2
        grad = barrier_grad_xcorr(v, R);
        %grad = -1./v;
    end
    if nargout >= 3
        hess = barrier_hess_xcorr(v, R);
        %hess = diag(1./(v.^2));
    end
end