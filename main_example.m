clear all; close all;
addpath('./Utils');
setGlobalBioImPath;
rng(1);

% Generate data points
[y, z] = generate_data(50, 5, 0.2);


%% Example for single regularization parameter

lambda = 7*1e-3;

% Compute parametric form of sparsest solution
[a_sol, x_sol, p_sol] = gBLASSO_sol(y, z, lambda);

% Function handle of sparsest solution
sparsest_sol = @(t) linear_spline(t, a_sol, x_sol, p_sol);

% Plot sparsest solution
font_size = 15; line_width = 2; marker_size = 12;
margin = (y(end) - y(1))/10; xmin = y(1) - margin; xmax = y(end) + margin;
t_grid = linspace(xmin, xmax, 1e4);
% Plot sparsest solution
figure;
plot(y, z, 'kx','LineWidth', line_width, 'Markersize', marker_size);
leg = {'Data points'}; ax = gca; set(ax, 'FontSize', font_size); hold on;
plot(t_grid, sparsest_sol(t_grid), 'Linewidth', line_width);
leg = [leg, {'Sparsest solution'}];

knots = x_sol(abs(a_sol) > 1e-5);
if ~isempty(x_sol)
    ax.ColorOrderIndex = 1; plot(x_sol, sparsest_sol(x_sol), 'o', 'LineWidth', ...
        line_width, 'Markersize', marker_size);
    leg = [leg, {'Knots'}];
end

if lambda ~= 0
    ax.ColorOrderIndex = 1;
    plot(y, sparsest_sol(y), 'x','LineWidth', line_width, 'Markersize', marker_size);
    leg = [leg, {'Denoised data points'}];
end
xlim([xmin, xmax]); legend(leg, 'Location', 'Best');

%% Example of UI with varying lambda

user_interface(y, z);
