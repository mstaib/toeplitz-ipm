function val = conjugate_barrier(s)

% feasible starting point
x = [4;4;2;1];
R = levinson_durbin(x);

for k=1:10
    val = -s'*x + barrier_value(x, R);

    grad = barrier_grad(x, R);
    hess = barrier_hess(x, R);

    full_grad = s + grad;

    if sum(sum(isinf(hess))) == 0 && sum(sum(isnan(hess))) == 0
        break
    end

    del_x = hess\full_grad;
    R_new = levinson_durbin(x - del_x);

    if sum(sum(isinf(R))) == 0 && sum(sum(isnan(R))) == 0
        x = x - del_x;
        R = R_new;
    else
        break
    end
end

val = -s'*x + barrier_value(x, R);

end
