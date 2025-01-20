import numpy as np
from matplotlib import pyplot as plt

from denoising import denoise_y
from sparsification import sparsest_interpolant, linear_spline

# Load data
data = np.loadtxt("../data.csv", delimiter=",")
x, y = data[:, 0], data[:, 1]

# Regularization parameter
lamb = 1e-2

# Compute denoised y
y_denoised = denoise_y(x, y, lamb, rho=lamb)

# Compute sparsest linear spline that connects denoised data points
knots, amplitudes, polynomial = sparsest_interpolant(x, y_denoised)

# Plot result
margin = (x[-1]-x[0]) / 10
t_grid = np.concatenate(([x[0]-margin], knots, [x[-1]+margin]))
fig = plt.figure()
ax = plt.gca()
ax.plot(x, y, 'x', label='Original data', markersize=10)
if lamb > 0:
    ax.plot(x, y_denoised, 'x', label='Denoised data', markersize=10)
ax.plot(t_grid, linear_spline(t_grid, knots, amplitudes, polynomial), label='Sparsest solution')
if len(knots) > 0:
    ax.plot(knots, linear_spline(knots, knots, amplitudes, polynomial), 'o', label='Knots')
ax.legend()
fig.show()
