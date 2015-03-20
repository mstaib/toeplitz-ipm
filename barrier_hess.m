function hess = barrier_hess(z, R)
n = length(z);

N = 2^ceil(log2(2*(n+1)));

hess_middle_sum = zeros(N);

% apparently .' is transpose and just ' is conjugate transpose

Rk_term_1 = zeros(N);
Rk_term_2 = zeros(N);
Rks = [];
for k=1:n
    Rk = fft(R(:,k), N);
    Rks = [Rks Rk];
    Rk_term_1 = Rk_term_1 + Rk * Rk.';
    Rk_term_2 = Rk_term_2 + Rk * Rk';
end

Rl_conj = conj(Rks(:,1));
lsum = diag(Rl_conj)*((Rk_term_1)*diag(Rl_conj) + Rk_term_2*diag(Rks(:,1)));
for l = 2:n
    Rl_conj = conj(Rks(:,l));
    lsum = lsum + diag(Rl_conj)*((Rk_term_1)*diag(Rl_conj) + Rk_term_2*diag(Rks(:,l)));
end
%Rl_term_1 = Rk_term_1.';
%Rl_term_2 = Rk_term_2';

%hess_middle_sum = Rk_term_1.*Rl_term_1 + Rk_term_2.*Rl_term_2;
lsum = lsum * fft(eye(N),N);

hess_tall = ifft(2*lsum, N);
hess = real(hess_tall(1:n,1:n));
end
