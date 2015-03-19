rand('state',1);

u = [4;4;2;1];

A = rand(3, 4);
b = A*u;

c = [2; 1; 1; 1];

[uopt, obj, uhist] = non_symm_primal_dual_ipm(c, A, b, u, 0.1, 0.5);
