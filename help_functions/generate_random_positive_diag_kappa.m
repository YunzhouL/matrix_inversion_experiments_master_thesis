% This function generates a random diagonal matrix D with
% k(D) = kappa with positive entrices

function Lambda = generate_random_positive_diag_kappa(n, kappa)
    % n: size of the matrix
    % kappa: desired condition number

    lambda = zeros(n, 1);

    % Step 1: set the extreme values with random signs
    lambda(1) = kappa;   % κ
    lambda(n) = 1;       % 1

    % Step 2: fill in λ₂ to λ_{n-1}
    for i = 2:n-1
        lambda(i) = 1 + (kappa - 1) * rand();   % uniform in [1, κ]
    end

    % Step 3: form the diagonal matrix
    Lambda = diag(lambda);
end