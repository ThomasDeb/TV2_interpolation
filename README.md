# TV2_interpolation
Regression method for 1D data that minimizes second-order total-variation regularization with the minimum number of piecewise-linear regions. This repository uses the [GlobalBioIm library](https://github.com/Biomedical-Imaging-Group/GlobalBioIm).

## Getting started

- gBLASSO_sol.m solves the (g-BLASSO) regression problem for a fixed regularization parameter
- user_interface.m does so with a GUI that allows to user to vary the regularization parameter manually
- main_example.m gives an example of how these functions are applied

## Reference

[Sparsest Piecewise-Linear Regression of One-Dimensional Data](https://www.sciencedirect.com/science/article/pii/S0377042721006130)  <br />
Journal of Computational and Applied Mathematics, vol. 406, p. 114044, May 2022. <br />
T. Debarre, Q. Denoyelle, M. Unser, and J. Fageot.

[Pocket Guide to Solve Inverse Problems with GlobalBioIm](https://iopscience.iop.org/article/10.1088/1361-6420/ab2ae9)  <br />
Inverse Problems, vol. 35, no. 10, pp. 10-20, October 2019. <br />
E. Soubies, F. Soulez, M. McCann,  T-a. Pham, L. Donati, T. Debarre, D. Sage, and M. Unser.
