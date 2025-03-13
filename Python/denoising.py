from typing import Tuple

import numpy as np
import scipy as sp


def denoise_y(x: np.ndarray, y: np.ndarray, lamb: float, rho: float = 1.0, max_iter: int = int(1e4),
                       relative_tol: float = 1e-3) -> np.ndarray:
    """ Solve Problem (36) in [1]_ using ADMM [2]_.

    Parameters
    ----------
    x : ndarray
        Array of x-coordinates of data points
    y : ndarray
        Array of y-coordinates of data points
    lamb : float
        Regularization parameter lambda
    rho : float
        Internal ADMM parameter [2]_
    max_iter : int
        Maximum number of iterations for ADMM
    relative_tol : float
        Tolerance parameter for ADMM stopping criterion (:math:`\epsilon^\mathrm{abs}` in [2]_)

    Returns
    -------
    y_lambda : ndarray
        Array of denoised y-coordinates of data points (:math:`\mathbf{y}_\lambda` in [1]_)

    References
    ----------
    .. [1] Debarre, T., Denoyelle, Q., Unser, M., and Fageot, J. "Sparsest Continuous Piecewise-Linear Spline
           Representation of One-Dimensional Data." Journal of Computational and Applied Mathematics, 2022.

    .. [2] Boyd, S., et al. “Distributed Optimization and Statistical Learning via the Alternating Direction
           Method of Multipliers.” Foundations and Trends in Machine Learning, 2011.
    """
    if x.size != y.size:
        raise Exception("x and y must be of the same size")
    if x.size < 2:
        raise Exception("Must have at least 2 data points")

    lamb_max, polynomial = lambda_max(x, y)
    if lamb >= lamb_max:
        # If lamb is too high, the problem amounts to linear regression
        y_denoised = polynomial[0] * np.ones_like(x) + polynomial[1] * x
    elif lamb > 0:
        # Otherwise, solve denoising problem using ADMM
        # Define matrices
        L = _regularization_matrix(x)  # Regularization matrix
        Lt = L.transpose()
        M = sp.sparse.identity(len(x)) + rho * (Lt @ L)  # Matrix to be inverted at every x-update step in ADMM
        M_diags = _diags_banded(M)  # Stacked diagonals of M to use Scipy's solveh_banded()
        # ADMM initialization
        xk, zk, yk = y, L @ y, np.zeros(L.shape[0])
        # Run ADMM
        for it in range(max_iter):
            z_prev = zk  # Needed for stopping criterion
            # ADMM updates (notations follow Boyd et al. 2011)
            b = y + rho * Lt @ (zk - yk / rho)
            xk = sp.linalg.solveh_banded(M_diags, b, lower=True)
            Lxk = L @ xk  # Precomputation
            zk = _prox_L1(Lxk + yk / rho, lamb / rho)
            yk += rho * (Lxk - zk)

            # ADMM stopping criterion (as suggested by Boyd et al. 2011)
            primal_res = np.linalg.norm(Lxk - zk)
            dual_res = rho * np.linalg.norm(Lt @ (zk - z_prev))
            primal_eps = relative_tol * max(np.linalg.norm(Lxk), np.linalg.norm(zk))
            dual_eps = relative_tol * np.linalg.norm(Lt @ yk)
            if (primal_res <= primal_eps) and (dual_res <= dual_eps):
                break
        y_denoised = xk

    else:
        # If lamb = 0, no denoising
        y_denoised = y
    return y_denoised


def lambda_max(x: np.ndarray, y: np.ndarray) -> Tuple[float, np.ndarray]:
    """ Compute maximum regularization parameter lambda above which the problem amounts to linear regression.
     Returns maximum lambda and the linear regression solution. """
    if x.size != y.size:
        raise Exception("x and y must be of the same size")
    n = len(x)

    # Compute parameters of polynomial (solution of a 2x2 system)
    det = n * np.sum(x ** 2) - np.sum(x) ** 2
    polynomial = (1 / det) * np.array([[np.sum(x ** 2), -np.sum(x)], [-np.sum(x), n]]).dot(
        np.array([np.sum(y), np.dot(x, y)]))

    if n > 2:
        h = y - (polynomial[0] * np.ones_like(x) + polynomial[1] * x)
        lamb_max = max(np.abs(x[1:-1] * np.cumsum(h)[:-2] - np.cumsum(h * x)[:-2]))
    elif n == 2:
        lamb_max = 0  # If 2 data points, lambda is irrelevant (the solution is always linear regression)
    return lamb_max, polynomial


def _regularization_matrix(x: np.ndarray) -> sp.sparse.diags:
    """ Compute the L matrix defined in Equation (37) in Debarre et al. 2022."""
    M = len(x)
    v = 1 / (x[1:] - x[:-1])
    return sp.sparse.diags([v[:-1], -(v[:-1] + v[1:]), v[1:]], [0, 1, 2], shape=(M-2, M))


def _diags_banded(M: sp.sparse.sparray) -> np.ndarray:
    """ Return stacked diagonals of banded M matrix to apply solveh_banded. """
    diag = M.diagonal()
    diag1 = M.diagonal(k=1)
    diag2 = M.diagonal(k=2)
    return np.stack((diag, np.concatenate([diag1, [0]]), np.concatenate([diag2, [0, 0]])))


def _prox_L1(x: np.ndarray, sigma: float):
    """ Compute proximal operator of the L1 norm. """
    return np.sign(x) * np.maximum(np.abs(x) - sigma, 0)
