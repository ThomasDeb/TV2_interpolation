function [y, z, gt] = generate_data(M, num_knots_GT, noise)
% Generate data points taken from a linear spline GT

a = randn(num_knots_GT, 1); a = (a+sign(a))/2;
x = (rand(num_knots_GT, 1) + (0:num_knots_GT-1)')/num_knots_GT;  p = [0, 0];
gt = @(t) linear_spline(t, a, x, p);
y = (rand(M, 1) + (0:M-1)')/M; z = gt(y) + noise*randn(M, 1);
end

