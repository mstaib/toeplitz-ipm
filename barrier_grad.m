function grad = barrier_grad(z, R)
n = length(z);

N = 2*n;

grad_sum = zeros(N, 1);
for k=1:n
    Rk = fft(R(:,k), N);
    grad_sum += Rk.*conj(Rk);
end

grad = -2/N * ifft(grad_sum, N);

end
