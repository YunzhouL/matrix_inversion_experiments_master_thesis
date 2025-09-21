% This function computes M^{-1} for a complex matrix M,
% with only real matrix operations
% M = A + iB then
% M^{-1} = (A+BA^{-1}B)^{-1} - iA^{-1}B * (A+BA^{-1}B)^{-1}

function X = inverse_Frobenius(M)
    A = real(M);
    B = imag(M);
    X_2 = A_inverse_B(A, B);            % Compute A^{-1}B
    X_3 = inverse_via_LU(A + B*X_2);    % Compute (A+BA^{-1}B)^{-1}
    X_4 = X_2*X_3;                      % Compute A^{-1}B * (A+BA^{-1}B)^{-1}
    X = X_3 - 1i*X_4;
end