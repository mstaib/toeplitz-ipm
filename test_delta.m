rand('state',1);

A = rand(3, 4);
b = A*[8;4;2;1];

c = ones(4,1);

u = [4;4;2;1];
R = levinson_durbin(u);
grad = barrier_grad(u, R);
hess = barrier_hess(u, R);

t = 0.5;
[d, y] = delta(c, A, b, grad, hess, t);
