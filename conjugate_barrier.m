function val = conjugate_barrier(s)

% feasible starting point
x = [4;4;2;1];

for k=1:10
    R = levinson_durbin(x);
    val = -s'*x + barrier_value(x, R);

    grad = barrier_grad(x, R);
    hess = barrier_hess(x, R);
    x
    s


    full_grad = s + grad;

    x = x - hess\full_grad;
end

R = levinson_durbin(x);
val = -s'*x + barrier_value(x, R);

end
