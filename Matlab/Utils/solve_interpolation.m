function [a_sol, x_sol, p_sol] = solve_interpolation(y, z, sparsity_tol)
% Outputs parametrizations for both solution and certificate of the interpolation
% problem s(y(i)) = z(i);

[y, I] = sort(y); z = z(I);
if size(y, 1) == 1
    y = y';
end
if size(z, 1) == 1
    z = z';
end

M = length(y);
x_sol = y(2:end-1);
[a_sol, p_sol] = connect_points(y, z);

x_sol_nz = x_sol(abs(a_sol) > sparsity_tol);
a_sol_nz = a_sol(abs(a_sol) > sparsity_tol);

% Find "phantom knots" (a=0) which are in saturation zones: they must be kept in
% mind, whereas the others can be discarded
x_sol_z = [];
a_sol_z = [];
signs = (abs(a_sol) > sparsity_tol) .* sign(a_sol);
nz_indices = find(signs);
if length(nz_indices) > 1
    idx = 1;
    for i = nz_indices(1) + 1 : nz_indices(end) - 1
        if (signs(i) == 0 && signs(nz_indices(idx)) * signs(nz_indices(idx+1)) == 1)
            % Keep phantom knots in saturation zones
            x_sol_z = [x_sol_z; x_sol(i)];
            a_sol_z = [a_sol_z; a_sol(i)];
        elseif (signs(i) ~= 0)
            idx = idx + 1;
        end
    end
end


K = length(x_sol); % Sparsity of canonical solution

% Compute relevant knots (including phantom knots in saturation zones)
[x_sol_z_nz, I] = sort([x_sol_z; x_sol_nz]);
a_sol_z_nz = [a_sol_z; a_sol_nz];
a_sol_z_nz = a_sol_z_nz(I);
a_sol_z_nz_tol = (abs(a_sol_z_nz) > sparsity_tol) .* a_sol_z_nz;
K = length(x_sol_z_nz);

a_sol_temp = [];
x_sol_temp = [];
k = 1;
while k <= K
    n = 0; % Number of consecutive saturations after knot x_k
    while ((k+n+1) <= K) && (a_sol_z_nz_tol(k+n+1) * a_sol_z_nz_tol(k) >= 0)
        n = n+1;
    end
    for i = 0 : ceil(n/2) - 1
        % Add barycenter knot
        a_sol_temp = [a_sol_temp; a_sol_z_nz(k+2*i) + a_sol_z_nz(k+2*i+1)];
        x_sol_temp = [x_sol_temp; (a_sol_z_nz(k+2*i) * x_sol_z_nz(k+2*i) + a_sol_z_nz(k+2*i+1) * x_sol_z_nz(k+2*i+1)) / a_sol_temp(end)];
    end
    if mod(n, 2) == 0
        % Keep existing knot
        a_sol_temp = [a_sol_temp; a_sol_z_nz(k+n)];
        x_sol_temp = [x_sol_temp; x_sol_z_nz(k+n)];
    end
    k = k+n+1;
end
a_sol_temp = sparsify_amplitudes(a_sol_temp, sparsity_tol);
x_sol = x_sol_temp(a_sol_temp ~= 0);
a_sol = a_sol_temp(a_sol_temp ~= 0);




