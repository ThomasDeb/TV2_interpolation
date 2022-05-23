clear all; close all;
setGlobalBioImPath;

[y, z, gt] = generate_data(50, 5, 0.2);

% Data points: (y(i), z(i))


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
sp_label = param.sp_label;
% ADMM parameters
rho = lambda/10; max_iter = 1000; relative_tol = 1e-14; iterVerb = 1; iterNum = 1;
tic;
[z_sol] = compute_z_sol(y, z, lambda, rho, max_iter, relative_tol, iterVerb, iterNum);
toc;

% Construct solution and certificate
sparsity_tol = 1e-5;
tic;
[a_sol, x_sol, p_sol, ~, ~] = solve_interpolation(y, z_sol, sparsity_tol);
toc;
knots = x_sol(abs(a_sol) > sparsity_tol);
sparsity = sum(abs(a_sol) > sparsity_tol);
solution = @(t) linear_spline(t, a_sol, x_sol, p_sol);

% Plot
font_size = 15; line_width = 2; marker_size = 12;
% Plot sparsest solution
plot(ax, y, z, 'kx','LineWidth', line_width, 'Markersize', marker_size);
leg = {'Original data points'}; set(ax, 'FontSize', font_size); hold(ax, 'on');
fplot(ax, solution, param.x_lim, 'Linewidth', line_width);
leg = [leg, {'Solution'}];
if ~isempty(knots)
    ax.ColorOrderIndex = 1; plot(ax, knots, solution(knots), 'o', 'LineWidth', ...
        line_width, 'Markersize', marker_size);
    leg = [leg, {'Knots'}];
end

if lambda ~= 0
    ax.ColorOrderIndex = 1;
    plot(ax, y, z_sol, 'x','LineWidth', line_width, 'Markersize', marker_size);
    leg = [leg, {'Measurements of solution set'}];
    sp_label.Text = sprintf('Lambda = %.1e \nSparsity: %i', lambda, sparsity);
else
    sp_label.Text = sprintf('Constrained problem \nSparsity: %i', sparsity);
end
legend(ax, leg,'Interpreter','latex', 'Location', 'Best');
ax.XLim = param.x_lim;
if lambda == 0
    ax.YLim = [min([param.y_lim(1), ax.YLim(1)]), max([param.y_lim(2), ax.YLim(2)])];
else
    ax.YLim = param.y_lim;
end
end



