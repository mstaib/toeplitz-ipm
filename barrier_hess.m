function hess = barrier_hess(z, Rk)
    n = length(z);
    N = size(Rk,1);
    % apparently .' is transpose and just ' is conjugate transpose
    Rk_term_1 = zeros(N);
    Rk_term_2 = zeros(N);
    for k=1:n
        Rk_term_1 = Rk_term_1 + Rk(:,k) * Rk(:,k).';
        Rk_term_2 = Rk_term_2 + Rk(:,k) * Rk(:,k)';
    end
    Rl_conj = conj(Rk(:,1));
    lsum = diag(Rl_conj)*((Rk_term_1)*diag(Rl_conj) + Rk_term_2*diag(Rk(:,1)));
    for l = 2:n
        Rl_conj = conj(Rk(:,l));
        lsum = lsum + diag(Rl_conj)*((Rk_term_1)*diag(Rl_conj) + Rk_term_2*diag(Rk(:,l)));
    end
    lsum = lsum * fft(eye(N),N);
    hess_tall = ifft((2/N)*lsum, N);
    hess = real(hess_tall(1:n,1:n));
end
