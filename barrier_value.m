function val = barrier_value(z, R)
% R comes from Levinson-Durbin, such that inv(toeplitz(z)) = R*R'

val = 2 * sum(log(diag(R)));
