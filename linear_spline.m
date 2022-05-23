function [y] = linear_spline(t, a, x, p)
% y = p(t) + sum a_k rho(t - x_k)

y = 0 * t;
if nargin > 3
    y = p(1) + p(2) * t;
end
for k = 1 : length(a)
y = y + a(k) * ((t - x(k)) > 0) .* (t - x(k));
end

end

