function hess = barrier_hess(z, R)
n = length(z);

N = 2*n;

hess_middle_sum = zeros(N);

% apparently .' is transpose and just ' is conjugate transpose

Rk_term_1 = zeros(N);
Rk_term_2 = zeros(N);
for k=1:n
    Rk = fft(R(:,k), N);

    Rk_term_1 = Rk_term_1 + Rk * Rk';
    Rk_term_2 = Rk_term_2 + Rk * Rk.';
end

Rl_term_1 = Rk_term_1.';
Rl_term_2 = Rk_term_2';

hess_middle_sum = Rk_term_1.*Rl_term_1 + Rk_term_2.*Rl_term_2;

hess_tall = real(2/N * ifft(hess_middle_sum * fft(eye(n), N), N));
hess = hess_tall(1:n, :);

end
