function grad = barrier_grad(z, R)
n = length(z);

N = 2*n;

grad_sum = zeros(N, 1);
for k=1:n
    Rk = fft(R(:,k), N);
    grad_sum = grad_sum + Rk.*conj(Rk);
end

grad_ifft = ifft(grad_sum, N);

grad = -2 * grad_ifft(1:n);

end
