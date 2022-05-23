clear all; close all;
setGlobalBioImPath;

% Data points are (y(i), z(i))
 y = [0; 1; 2; 3; 4];
 z = [0; 2; 3; 3; 2];

 % Regularization parameter (lambda = 0 -> constrained interpolation problem)
lambda = 0; 

% ADMM parameters
rho = 1e0; max_iter = 1000; relative_tol = 1e-14; iterVerb = 1; iterNum = 1;

% Compute "denoised" data points
[z_sol] = compute_z_sol(y, z, lambda, rho, max_iter, relative_tol, iterVerb, iterNum);

%% Construct solution and certificate
sparsity_tol = relative_tol;

[a_sol, x_sol, p_sol, a_cert, x_cert] = solve_interpolation(y, z_sol, sparsity_tol);
knots = x_sol(abs(a_sol) > sparsity_tol);
sparsity = sum(abs(a_sol) > sparsity_tol);
solution = @(t) linear_spline(t, a_sol, x_sol, p_sol);

fprintf('Minimum sparsity (lambda = %e): %i\n', lambda, sparsity);

%% Plots
font_size = 15; line_width = 2; marker_size = 12;
margin = (y(end) - y(1))/10; xmin = y(1) - margin; xmax = y(end) + margin; 
% Plot sparsest solution
figure(1);
plot(y, z, 'kx','LineWidth', line_width, 'Markersize', marker_size);
leg = {'Original data points'}; ax = gca; set(ax, 'FontSize', font_size); hold on;
fplot(solution, [xmin, xmax], 'Linewidth', line_width);
leg = [leg, {'Solution'}];
if ~isempty(knots)
    ax.ColorOrderIndex = 1; plot(knots, solution(knots), 'o', 'LineWidth', ...
        line_width, 'Markersize', marker_size);
    leg = [leg, {'Knots'}];
end

if lambda ~= 0
    ax.ColorOrderIndex = 1;
    plot(y, z_sol, 'x','LineWidth', line_width, 'Markersize', marker_size);
    leg = [leg, {'Measurement of solution set'}];
end



