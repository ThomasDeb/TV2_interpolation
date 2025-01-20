function [a_sol, x_sol, p_sol] = gBLASSO_sol(y, z, lambda, sparsity_tol, rho, max_iter, relative_tol, iterVerb, iterNum)
% Output sparsest solution of (g-BLASSO)
%
% Inputs:
% - y, z        : vectors of x- and y-coordinates of data points
% - lambda      : regularization parameter
% - sparsity_tol: if the magnitude of a weight is below this threshold, the knot is deleted
% - rho, max_iter, relative_tol, iterVerb, iterNum: ADMM parameters (see GlobalBioIm documentation)
%
% Outputs:
% - a_sol: vector of weights of the sparsest solution
% - x_sol: vector of knots of the sparsest solution
% - p_sol: polynomial term of the sparsest solution

% ADMM parameters
if nargin < 9
    iterNum = 1;
    if nargin < 8
        iterVerb = 100;
        if nargin < 7
            relative_tol = 1e-14;
            if nargin < 6
                max_iter = 2000;
                if nargin < 5
                    rho = 1e0;
                    if nargin < 4
                        sparsity_tol = 1e-5;
                    end
                end
            end
        end
    end
end

% Compute denoised data points
z_sol = compute_z_sol(y, z, lambda, rho, max_iter, relative_tol, iterVerb, iterNum);

% Compute sparsest solution
[a_sol, x_sol, p_sol] = solve_interpolation(y, z_sol, sparsity_tol);

end