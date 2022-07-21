function [a, p] = connect_points(y, z)
% Outputs weights and null space element of D^ splines which connects
% points (y_i, z_i). No knots at y(1) and y(M).

alpha = (z(2:end) - z(1:end-1)) ./ (y(2:end) - y(1:end-1)); % List of slopes
p = [z(1) - alpha(1) * y(1); alpha(1)]; % Null space coefficients
a = alpha(2:end) - alpha(1:end-1); % Weights


end

