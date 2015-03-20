function grad = barrier_grad(z, Rk)
    n = length(z);
    N = size(Rk,1);
    grad_sum = -2*sum(abs(Rk).^2,2);
    grad = ifft(grad_sum, N);
    grad = real(grad(1:n));
end
