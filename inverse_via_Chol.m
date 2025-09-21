% This function computes for hpd or spd matrix M its inverse M^{-1}
% using Cholesky decomposition and backsubstitutions, the flop count is 2n^3

function X = inverse_via_Chol(M)
    L = chol(M, "lower");
    X1 = inv(L');    % Solve X1 = L^{-1}
                    % The alternative way computing X1 is not ideal:
                    % X1 = linsolve(U, eye(size(M), struct('UT', true));
                    
    X2 = linsolve(L', X1', struct('UT', true)); % Solve X2 = (X1 * L^{-1})' 
    X = X2';
end