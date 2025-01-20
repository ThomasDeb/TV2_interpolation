function [L] = regularization_matrix(y)
% L matrix for optimization problem min ||y - z||_2^2 + \lambda ||L z||_1
M = length(y);
L = zeros(M-2, M);
inv_dy = 1 ./ (y(2:end) - y(1:end-1));

for i = 1 : M-2
    L(i, i) = inv_dy(i);
    L(i, i+1) = - (inv_dy(i+1) + inv_dy(i));
    L(i, i+2) = inv_dy(i+1);
end
end


