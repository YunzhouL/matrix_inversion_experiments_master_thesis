% This function computes M^{-1} for a HPD complex matrix M,
% with only real matrix operations
% M = A + iB then
% M^{-1} = (A+BA^{-1}B)^{-1} - iA^{-1}B * (A+BA^{-1}B)^{-1}

function X = inverse_Frobenius_Chol(M)
    A = real(M);
    B = imag(M);
    L1 = chol(A, "lower");
    X1 = linsolve(L1, B, struct('LT', true));       % X1 = L1^{-1}B
    X2 = linsolve(L1', X1, struct('UT', true));     % X2 = A^{-1}B
    X3 = inverse_via_Chol(A - X1' * X1);
    X4 = X2*X3;                      % Compute A^{-1}B * (A+BA^{-1}B)^{-1}
    X = X3 - 1i*X4;
end