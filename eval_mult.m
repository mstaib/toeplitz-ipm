function vals = eval_mult(func, points)

vals = zeros(size(points));
for kk=1:length(points)
    vals(kk) = func(points(kk));
end

end
