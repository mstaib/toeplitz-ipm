function R = levinson_durbin(z)
% produces lower triangular R so that inv(toeplitz2(z)) = R*R'

z = z(:);
n = length(z);
E = diag(ones(n-1,1), -1);

R = zeros(n);

R(1,1) = 1/sqrt(2*z(1)); %1/sqrt(2*z(1)); % make this 2*z(1) if you want diagonal to be 2*z(1)

for k=1:(n-1)
    rk_tilde = [flipud(R(1:k,k)); zeros(n-k,1)];
    alpha_k = -R(k,k) * z' * E * R(:,k);

    R(:,k+1) = 1/sqrt(1-alpha_k^2) * (E*R(:,k) + alpha_k*rk_tilde);
end

end
