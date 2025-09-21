% This function is used in conjunction with generate_two_random_diag_kappa_with_relation 
% to create test matrices while maintaining constant κ(A), κ(B), and κ(M).
% The two functions work together to generate matrices with either 
% large κ(A+BA⁻¹B) or small κ(A+BA⁻¹B)

function [A, B] = generate_comp_greater(n, kappa)
    % Generate two diagonal matrices A and B
    % Both have condition number = kappa and κ(A+BA^{-1}B) = (kappa^2+1)/2

    a_vals = zeros(n, 1);
    b_vals = zeros(n, 1);

    % Step 1: set extreme values with random signs
    a_vals(1) = kappa * (-1)^(randi([0, 1])); 
    a_vals(n) = 1 * (-1)^(randi([0, 1]));
    a_vals(2) = 1 * (-1)^(randi([0, 1])); 

    b_vals(1) = kappa * (-1)^(randi([0, 1])); 
    b_vals(n) = 1 * (-1)^(randi([0, 1]));
    b_vals(2) = kappa * (-1)^(randi([0, 1])); 

    % Step 2: fill λ₂ to λ_{n-1} with the constraint
    for j = 3:n-1
        valid = false;
        while ~valid
            % Random Aj in ±[1, κ]
            a = 1 + (kappa - 1) * rand();
            if rand() < 0.5
                a = -a;
            end

            % Random Bj in ±[1, κ]
            b = 1 + (kappa - 1) * rand();
            if rand() < 0.5
                b = -b;
            end

            % Check constraint
            val = a + (b^2) / a;
            if (abs(val) > 2 && abs(val) < 2 * kappa)
                valid = true;
                a_vals(j) = a;
                b_vals(j) = b;
            end
        end
    end

    % Step 3: build diagonal matrices
    A = diag(a_vals);
    B = diag(b_vals);
end