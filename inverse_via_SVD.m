% This function computes M^{-1} using SVD, the flop count
% is approximately 19n^3

function X = inverse_via_SVD(M)
    [U, D, V] = svd(M);
    X = V * diag(1 ./ diag(D)) * U';    % Solve X = V * D^{-1} * U'
end