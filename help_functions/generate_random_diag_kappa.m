% This function generates a random diagonal matrix D with κ(D) = kappa

function Lambda = generate_random_diag_kappa(n, kappa)
    % n: size of the matrix
    % kappa: desired condition number

    lambda = zeros(n, 1);

    % Step 1: set the extreme values with random signs
    lambda(1) = kappa * (-1)^(randi([0,1]));   % ±κ
    lambda(n) = 1 * (-1)^(randi([0,1]));       % ±1

    % Step 2: fill in λ₂ to λ_{n-1}
    for i = 2:n-1
        val = 1 + (kappa - 1) * rand();   % uniform in [1, κ]
        if rand() < 0.5
            val = -val;                  % flip sign with 50% probability
        end
        lambda(i) = val;
    end

    % Step 3: form the diagonal matrix
    Lambda = diag(lambda);
end