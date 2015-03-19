function val = conjugate_barrier(s)

val = -log(det(toeplitz2(s))) - length(s);

end
