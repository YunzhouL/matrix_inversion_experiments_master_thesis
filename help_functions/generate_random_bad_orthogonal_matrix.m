% This function generates Q, D, U corresponding to Example 4.10

function [Q,D,U] = generate_random_bad_orthogonal_matrix(n, delta)
    
    % Step 1: Generate random eigenvalues with |λ| = 1
    eigenvalues = [];

    theta = pi/2 + delta;
    a = cos(theta);
    b = sin(theta);
    eigenvalues = [eigenvalues, a + 1i*b, a - 1i*b];

    i = 3;
    while i <= n
        if i < n && rand > 0.95 % Randomly decide to add complex pair
            theta = 2 * pi * rand; % Random angle for e^{iθ}
            a = cos(theta);
            b = sin(theta);
            eigenvalues = [eigenvalues, a + 1i*b, a - 1i*b]; % Add a+ib, a-ib
            i = i + 2;
        else % Add real eigenvalue ±1
            eigenvalues = [eigenvalues, (-1)^randi([0, 1])];
            i = i + 1;
        end
    end
    
    % Step 2: Construct complex diagonal matrix D
    D = diag(eigenvalues);
    
    % Step 3: Generate random real orthogonal matrix P as basis
    A = randn(n, n);
    [P, ~] = qr(A); % P is real orthogonal, columns are real basis vectors
    
    % Construct complex unitary U with conjugate eigenvector pairs
    U = zeros(n, n) + 1i * zeros(n, n); % Initialize complex U
    pos_real = 1; % Position in real basis P
    pos_eig = 1;  % Position in eigenvalues and U
    while pos_eig <= n
        eig = eigenvalues(pos_eig);
        if isreal(eig) % Real eigenvalue (±1)
            U(:, pos_eig) = P(:, pos_real); % Real unit vector
            pos_real = pos_real + 1;
            pos_eig = pos_eig + 1;
        else % Complex conjugate pair
            u = P(:, pos_real);
            w = P(:, pos_real + 1);
            v = (u + 1i * w) / sqrt(2);     % Normalized eigenvector for a+ib
            vbar = (u - 1i * w) / sqrt(2);  % Conjugate for a-ib
            U(:, pos_eig) = v;
            U(:, pos_eig + 1) = vbar;
            pos_real = pos_real + 2;
            pos_eig = pos_eig + 2;
        end
    end
    
    % Step 4: Compute Q = U * D * U' (U' is Hermitian transpose)
    Q = U * D * U';
    
    % Step 5: Ensure Q is real (handle numerical precision)
    Q = real(Q);
    
    % Step 6: Verify orthogonality and reality
    if norm(Q' * Q - eye(n), 'fro') > 1e-10
        warning('Orthogonality check failed: Q''*Q is not identity');
    end
    if norm(imag(Q), 'fro') > 1e-10
        warning('Q is not real: imaginary part detected');
    end
end