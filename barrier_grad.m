function grad = barrier_grad(z, R)
n = length(z);
% smallest power of 2 greater than 2(n+1)
N = 2^ceil(log2(2*(n+1)));

grad_sum = zeros(N, 1);
for k=1:n
    Rk = fft(R(:,k), N);
    grad_sum = grad_sum + Rk.*conj(Rk);
end

grad = ifft(-2*grad_sum, N);
grad = grad(1:n);
grad(1) = grad(1)/2;
end
