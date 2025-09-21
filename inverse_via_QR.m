% This function computes M^{-1} using QR decomposition and
% a backsubstitution, the flop count is 11/3n^3

function X = inverse_via_QR(M)
    [Q, R] = qr(M);
    X = linsolve(R, Q', struct('UT', true));    % Solve X = R^{-1} * Q'
end