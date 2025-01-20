# TV2_interpolation
Matlab and Python implementation of a regression method for noisy 1D data that minimizes second-order total-variation regularization with the minimum number of piecewise-linear regions. The Matlab implementation uses the [GlobalBioIm library](https://github.com/Biomedical-Imaging-Group/GlobalBioIm) as a dependency. The Python version uses a custom implementation of the Alternating Direction Method of Multipliers (ADMM) algorithm; its only dependencies are Numpy and Scipy.

## Getting started

**Matlab**:
- gBLASSO_sol.m solves the (g-BLASSO) regression problem for a fixed regularization parameter
- user_interface.m does so with a GUI that allows to user to vary the regularization parameter manually
- main_example.m gives an example of how these functions are applied

**Python**:
- `denoise_y()` in denoising.py computes the denoised data points (with adjustable noise level) using ADMM
- `sparsest_interpolant()` in sparsification.py computes the linear spline that connects the noisy data points with the fewest possible knots using Algorithm 1 in Debarre et al. 2022
- main_example.py gives an example of how these functions are applied

## References

[Sparsest Piecewise-Linear Regression of One-Dimensional Data](https://www.sciencedirect.com/science/article/pii/S0377042721006130)  <br />
Journal of Computational and Applied Mathematics, vol. 406, p. 114044, May 2022. <br />
T. Debarre, Q. Denoyelle, M. Unser, and J. Fageot.

[Pocket Guide to Solve Inverse Problems with GlobalBioIm](https://iopscience.iop.org/article/10.1088/1361-6420/ab2ae9)  <br />
Inverse Problems, vol. 35, no. 10, pp. 10-20, October 2019. <br />
E. Soubies, F. Soulez, M. McCann,  T-a. Pham, L. Donati, T. Debarre, D. Sage, and M. Unser.

[Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers.](https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf)  <br />
Foundations and Trends® in Machine learning, vol. 3, no. 1, p. 1-122, 2011. <br />
S. Boyd, N. Parikh, E. Chu, B. Peleato, and J. Eckstein.
