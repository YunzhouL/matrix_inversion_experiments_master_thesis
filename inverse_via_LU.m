% This function computes M^{-1} using LU decomposition and 
% backsubstitutions, the flop count is 2n^3

function X = inverse_via_LU(M)
    [L, U, p] = lu(M,"vector");
    X1 = inv(U);    % Solve X1 = U^{-1}
                    % The alternative way computing X1 is not ideal:
                    % X1 = linsolve(U, eye(size(M), struct('UT', true));
                    
    X2 = linsolve(L', X1', struct('UT', true)); % Solve X2 = (X1 * L^{-1})' 
    p_inv1(p) = 1:numel(p);
    X2 = X2(p_inv1, :);
    X = X2';
end