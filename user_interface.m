function [] = user_interface(y, z, sparsity_tol, max_iter, relative_tol, verbose)
% User interface to plot the sparsest solution of the (g-BLASSO) with a
% varying regularization parameter
%
% Inputs:
% - y, z: vectors of x- and y-coordinates of data points
% - sparsity_tol: if the magnitude of a weight is below this threshold, the knot is not counted
% - rho, max_iter, relative_tol, verbose: ADMM parameters (see GlobalBioIm documentation)

if nargin < 3
    sparsity_tol = 1e-5;
    if nargin < 4
        max_iter = 2000;
        if nargin < 5
            relative_tol = 1e-14;
            if nargin < 6
                verbose = 0;
            end
        end
    end
end
param.sparsity_tol = sparsity_tol;
param.max_iter = max_iter; param.relative_tol = relative_tol;
param.iterVerb = 0; param.iterNum = 0; param.verbose = verbose;

margin = (y(end) - y(1))/10; xmin = y(1) - margin; xmax = y(end) + margin;
x_lim = [xmin, xmax]; param.x_lim = x_lim; param.y_lim = [0, 1];

[lamb_max, alpha] = lambda_max(y, z);
max_pow = ceil(log10(lamb_max)); min_pow = max_pow - 5;
tick_labels = strrep('1e_', '_', cellstr(strsplit(num2str(min_pow:max_pow))));

fig = uifigure;
ax = uiaxes(fig, 'Position', [20 80 500 320]);
sp_label = uilabel(fig, 'Position', [40 15 140 70], 'Fontsize', 14);
param.sp_label = sp_label; plot_solution(ax, 0, y, z, param); param.y_lim = ax.YLim;
checkbox = uicheckbox(fig, 'Value', 1, 'Text', 'Constrained problem', 'Position', [200 35 130 22], ...
    'ValueChangedFcn',@(checkbox, event) cBoxChanged(event, ax, y, z, param));
slider = uislider(fig, 'Position', [370 50 150 3], 'Limits', [min_pow, max_pow], 'ValueChangingFcn',...
    @(slider, event) sliderMoving(event, checkbox, ax, y, z, param), ...
    'MajorTicks', min_pow:max_pow, 'MajorTickLabels', tick_labels);
lambda_label = uilabel(fig, 'Position', [420 60 100 15], 'Text', 'Lambda');


% Compute axis limits
y_p = sort(alpha(1) + alpha(2)*x_lim);
y_lim = ax.YLim; y_lim = [min([y_lim(1), y_p(1)]), max([y_lim(2), y_p(2)])];
ax.XLim = x_lim; ax.YLim = y_lim;
param.y_lim = y_lim;
end

function  cBoxChanged(event, ax, y, z, param)
% if event.Value == 0
%     lambda = slider.Value;
% else
%     lambda = 0;
% end
hold(ax, 'off');
plot_solution(ax, 0, y, z, param);
end

function sliderMoving(event, checkbox, ax, y, z, param)

if checkbox.Value == 0
    lambda = 10^event.Value;
    hold(ax, 'off');
    plot_solution(ax, lambda, y, z, param);
end

end


function plot_solution(ax, lambda, y, z, param)
sp_label = param.sp_label; sparsity_tol = param.sparsity_tol;
% ADMM parameters
rho = lambda/10; % Rho adapted to lambda
[z_sol] = compute_z_sol(y, z, lambda, rho, param.max_iter, param.relative_tol, param.iterVerb, param.iterNum, param.verbose);

% Construct solution
[a_sol, x_sol, p_sol] = solve_interpolation(y, z_sol, sparsity_tol);
knots = x_sol(abs(a_sol) > sparsity_tol);
sparsity = sum(abs(a_sol) > sparsity_tol);
solution = @(t) linear_spline(t, a_sol, x_sol, p_sol);

% Plot
font_size = 15; line_width = 2; marker_size = 12;
% Plot sparsest solution
plot(ax, y, z, 'kx','LineWidth', line_width, 'Markersize', marker_size);
leg = {'Data points'}; set(ax, 'FontSize', font_size); hold(ax, 'on');
fplot(ax, solution, param.x_lim, 'Linewidth', line_width);
leg = [leg, {'Sparsest solution'}];
if ~isempty(knots)
    ax.ColorOrderIndex = 1; plot(ax, knots, solution(knots), 'o', 'LineWidth', ...
        line_width, 'Markersize', marker_size);
    leg = [leg, {'Knots'}];
end

if lambda ~= 0
    ax.ColorOrderIndex = 1;
    plot(ax, y, z_sol, 'x','LineWidth', line_width, 'Markersize', marker_size);
    leg = [leg, {'Denoised data points'}];
    sp_label.Text = sprintf('Lambda = %.1e \nSparsity: %i', lambda, sparsity);
else
    sp_label.Text = sprintf('Constrained problem \nSparsity: %i', sparsity);
end
legend(ax, leg, 'Location', 'Best');
ax.XLim = param.x_lim;
if lambda == 0
    ax.YLim = [min([param.y_lim(1), ax.YLim(1)]), max([param.y_lim(2), ax.YLim(2)])];
else
    ax.YLim = param.y_lim;
end
end
