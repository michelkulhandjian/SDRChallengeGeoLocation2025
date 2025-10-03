% Sparse Processing
function [STFT_sparse]  = SparseP (curOut, num_sparse_elem, totalFreqBins)
[M, I] = maxk(abs(curOut),num_sparse_elem,2);
x_ind = reshape(I', [],1);
z_val = reshape(M', [],1);
freq_ind = kron([1:totalFreqBins]', ones(1, num_sparse_elem));
y_ind = reshape(freq_ind', [],1);

Li = y_ind + (x_ind-1) * max(y_ind);
STFT_sparse = zeros(size(curOut));
STFT_sparse(Li) = curOut(Li);
end % Sparse P