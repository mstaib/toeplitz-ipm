function [ value, grad, hess ] = barrier(v)
    R = levinson_durbin(v);
    if nargout >=0
        value = barrier_value(v,R);
        %value = -sum(log(v));
    end
    if nargout >= 2
        % calculate DFT of cols of R
        n = length(v);
        N = 2^ceil(log2(2*n));
        Rk = zeros(N,n);
        for k=1:n
           Rk(:,k) = fft(R(:,k), N);
        end
        grad = barrier_grad(v, Rk);
        %grad2 = barrier_grad_xcorr(v,R);
        %assert(norm(grad-grad2) <= 1e-6);
        %grad = -1./v;
        if nargout >= 3
            hess = barrier_hess(v, Rk);
            %hess2 = barrier_hess_xcorr(v,R);
            %assert(norm(hess(:) - hess2(:)) <= 1e-6);
            %hess = diag(1./(v.^2));
        end
    end
end

