function [z_sol, ADMM] = compute_z_sol(y, z, lambda, rho, max_iter, relative_tol, iterVerb, iterNum, verbose)

if nargin < 9
    verbose = 1;
end

% Compute measurements of solution set
assert(issorted(y), 'The x-axis values must be sorted');
[lamb_max, alpha] = lambda_max(y, z);
ADMM = [];
if lambda == 0
    z_sol = z;
elseif lambda >= lamb_max
    p = @(x) alpha(1) + alpha(2) * x;
    z_sol = p(y);
else
    %% ADMM
    M = length(y);
    L = regularization_matrix(y);
    L_op = LinOpMatrix(L);
    C_fidelity = CostL2(M, z);
    C_regul = lambda * CostL1(L_op.sizeout);
    Cn = {C_regul};
    Hn = {L_op};
    solver = @(zz, rho, x) (eye(M)+rho*(L'*L)) \ (z + rho*L'*zz{1});
    ADMM = OptiADMM(C_fidelity, Cn, Hn, rho, solver);
    ADMM.OutOp = OutputOpti(true, iterVerb);
    ADMM.rho_n = rho;
    CvOp = TestCvgADMM(0, relative_tol);
    ADMM.CvOp = CvOp;
    ADMM.maxiter = max_iter;
    ADMM.CvOp.eps_rel = relative_tol;
    ADMM.ItUpOut = iterNum;
    ADMM.verbose = verbose;
    ADMM.run(y);
    z_sol = ADMM.xopt;
end


end

