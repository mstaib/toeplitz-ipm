function alpha = get_stepsize(x, s, y, del_x, del_s, del_y, beta, gamma)

n = length(x);

function val = function_to_zero(alpha)
    val = omega(x - alpha*del_x, s - alpha*del_s, y - alpha*del_y);
    val = val - omega(x, s, y);
    val = val - beta^2*gamma -(-gamma - log(1-gamma));
end

alpha_min = 0;
alpha_max = (1-beta)/(beta+sqrt(n));

top_end_low = 0;
top_end_high = alpha_max;

while imag(function_to_zero(top_end_high)) == 0
    top_end_high = top_end_high * 2;
end

while top_end_high - top_end_low >= 1e-6 && function_to_zero(top_end_low) <= 0
    mid = (top_end_high + top_end_low) / 2;
    val = function_to_zero(mid);
    if imag(val) ~= 0
        top_end_high = mid;
    else
        top_end_low = mid;
    end
end
%printf('exiting, low: %f, high: %f\n', top_end_low, top_end_high);

%points = 0:0.1:top_end_high;
%figure;
%plot(points, eval_mult(@(x) function_to_zero(x), points));

if function_to_zero(top_end_low) <= 0
    alpha = top_end_low;
else
    alpha_min = 0;
    alpha_max = top_end_low;
    while alpha_max - alpha_min >= 1e-6
        mid = (alpha_max + alpha_min)/2;
        if function_to_zero(mid) > 0
            alpha_max = mid;
        else
            alpha_min = mid;
        end
    end

    alpha = alpha_min;
end

if abs(function_to_zero(0)) < abs(function_to_zero(alpha))
    alpha = 0;
end

end
