function hess = barrier_hess_xcorr(z, R)
n = length(z);

hess = zeros(n);
for k=1:n
    for l=1:n
        xcorr_kl_long = xcorr(R(:,k), R(:,l));
        xcorr_lk_long = xcorr(R(:,l), R(:,k));

        xcorr_kl = xcorr_kl_long(n:end);
        xcorr_lk = xcorr_lk_long(n:end);

        hess = hess + xcorr_kl * (xcorr_lk + xcorr_kl)';
    end
end

hess = 2*hess;

end
