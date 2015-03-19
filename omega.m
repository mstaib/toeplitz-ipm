function val = omega(x, s, y)

n = length(x);
val = -log(det(toeplitz2(x))) + conjugate_barrier(s) + n*log((s'*x)/n) + n;

end
