function [amplitudes_sparse] = sparsify_amplitudes(amplitudes, sparsity_tol)
% Set amplitudes below tolerance (phantom knots) to exactly zero and compensate
% so that the solution has the exact same expression outisde areas with phantom knots

amplitudes_sparse = amplitudes;
zero_idx = find(abs(amplitudes) <= sparsity_tol);
amplitudes_sparse(zero_idx) = 0;
i = 1;
while i <= length(zero_idx)
    amplitudes_sparse(zero_idx(i)-1) = amplitudes_sparse(zero_idx(i)-1) + amplitudes_sparse(zero_idx(i));
    if i == length(zero_idx)
        break;
    end
    j = 0;
    while (i+j+1 <= length(zero_idx)) && (zero_idx(i+j+1) == zero_idx(i+j) + 1)
        amplitudes_sparse(zero_idx(i)-1) = amplitudes_sparse(zero_idx(i)-1) + amplitudes_sparse(zero_idx(i+j+1));
        j = j + 1;
    end
    i = i + j + 1;
end