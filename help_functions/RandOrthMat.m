function O = RandOrthMat(n)

% Generate random orthogonal matrix

    X = randn(n);
    [Q,R] = qr(X);
    R = diag(diag(R)./abs(diag(R)));
    O = Q*R;

end