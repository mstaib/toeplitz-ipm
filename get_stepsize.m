function alpha = get_stepsize(x, s, y, del_x, del_s, del_y, beta, gamma)

n = length(x);

function val = function_to_zero(alpha)
    val = omega(x - alpha*del_x, s - alpha*del_s, y - alpha*del_y);
    val = val - omega(x, s, y);
    val = val - beta^2*gamma -(-gamma - log(1-gamma));
end

alpha_min = 0;
alpha_max = (1-beta)/(beta+sqrt(n));

function_to_zero(0);
function_to_zero((1-beta)/(beta+sqrt(n)));

% find this by bisection
alpha = 0.01;

end
