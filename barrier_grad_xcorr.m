function grad = barrier_grad_xcorr(z, R)
n = length(z);

grad = zeros(n,1);
for k=1:n
    xcorr_k = xcorr(R(:,k), R(:,k));
    grad = grad + xcorr_k(n:end);
end

grad = -2 * grad;

end
