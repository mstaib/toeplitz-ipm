function mat = toeplitz2(z)

z_col = z(:);
z_col(1) = 2*z_col(1);
mat = toeplitz(z_col);
