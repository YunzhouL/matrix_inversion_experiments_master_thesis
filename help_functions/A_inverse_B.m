% This function computes A^{-1}B using LU decomposition and
% backsubstitutions, the flop count is 8/3n^3

function X = A_inverse_B(A, B)
    [L, U, p] = lu(A,"vector");
    B1 = B(p, :);   % Compute B1 = P * B
    X1 = linsolve(L, B1, struct('LT', true));   % Solve L * X1 = B1
    X = linsolve(U, X1, struct('UT', true));    % Solve U * X = X1
end