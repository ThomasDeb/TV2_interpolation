function [lambda, alpha] = lambda_max(y, z)
% Compute limit value of lambda above which the solution is a linear
% function
M = length(y);

% Compute this linear solution
A = [M, sum(y); ...
    sum(y), sum(y.^2)]; %<\nu(p_i), \nu(p_j)>
alpha = A \ [sum(z); sum(z .* y)]; % p = \sum alpha_i p_i
% p is the only element of the null space such that z - \nu(p) \in
% \nu(N_L) orthogonal

h = z - (alpha(1) * ones(M, 1) + alpha(2) * y); %z - \nu(p)
temp1 = cumsum(h);
temp2 = cumsum(h .* y);

lambda = max(abs(y(2:end-1) .* temp1(1:end-2) - temp2(1:end-2))); % ||\nu^*(temp)||_\infty



end

