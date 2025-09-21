% This program extends Test_5_2_2_part_1, demonstrating that using
% left/right relative residual for error evaluation generally yields
% results consistent with relative forward error measurements

% Get the current project directory and add helper functions to the path
projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'help_functions'));

% Define matrix sizes to test (from 100 to 2000 in steps of 100)
matrix_sizes = 100:100:2000;  % 20 different sizes in total
num_trials = 10;              % Number of trials for each matrix size
num_methods = 4;              % Number of methods to compare
average_errors = zeros(length(matrix_sizes), num_methods); % Store average errors

anomalies = {};  % Initialize cell array to record anomalies: method, matrix, time, average time, deviation, condition number

for idx = 1:length(matrix_sizes)
    n = matrix_sizes(idx);
    errors = zeros(num_trials, num_methods);

    % Perform multiple trials for statistical reliability
    for trial = 1:num_trials
        M = eye(n);
        rng(trial);  % Set fixed random seed for reproducibility

        % Generate M = A +iB with κ(M) = κ(A) = κ(B) = κ(A + BA⁻¹B) = 10
        while cond(M) <= 1 || cond(M) >= 100
            U1 = RandOrthMat(n);
            V1 = RandOrthMat(n);
            [S1, S2] = generate_two_random_diag_kappa_with_relation(n, 10);
            A = U1 * S1 * V1';
            B = U1 * S2 * V1';

            M = A + 1i * B;
        end

        norm_M = norm(M, 'fro');

        % Left relative residuals
        % Algo-LU
        X1 = inverse_via_LU(M);
        errors(trial, 1) = norm(X1 * M - eye(n), 'fro') / (norm(X1, 'fro') * norm_M);

        % Algo-Frob
        X2 = inverse_Frobenius(M);
        errors(trial, 2) = norm(X2 * M - eye(n), 'fro') / (norm(X2, 'fro') * norm_M);

        % Algo-QR
        X3 = inverse_via_QR(M);
        errors(trial, 3) = norm(X3 * M - eye(n), 'fro') / (norm(X3, 'fro') * norm_M);

        % Algo-SVD
        X4 = inverse_via_SVD(M);
        errors(trial, 4) = norm(X4 * M - eye(n), 'fro') / (norm(X4, 'fro') * norm_M);

        % Right relative residuals
        % X1 = inverse_via_LU(M);
        % errors(trial, 1) = norm(M * X1 - eye(n), 'fro') / (norm(X1, 'fro') * norm_M);
        % 
        % X2 = inverse_Frobenius(M);
        % errors(trial, 2) = norm(M * X2 - eye(n), 'fro') / (norm(X2, 'fro') * norm_M);
        % 
        % X3 = inverse_via_QR(M);
        % errors(trial, 3) = norm(M * X3 - eye(n), 'fro') / (norm(X3, 'fro') * norm_M);
        % 
        % X4 = inverse_via_SVD(M);
        % errors(trial, 4) = norm(M * X4 - eye(n), 'fro') / (norm(X4, 'fro') * norm_M);
    end

    % Calculate average errors for each method
    avg_error = mean(errors, 1);
    average_errors(idx, :) = avg_error;

    % Check for anomalies (trials where error exceeds 900% of average)
    for method = 1:num_methods
        for trial = 1:num_trials
            e = errors(trial, method);
            if e > (num_trials-1) * avg_error(method)
                condM = cond(M);
                anomalies{end+1} = struct( ...
                    'matrix_size', n, ...
                    'method', method, ...
                    'trial_error', e, ...
                    'avg_error', avg_error(method), ...
                    'percent_over', 100 * (e - avg_error(method)) / avg_error(method), ...
                    'cond_num', condM ...
                );
            end
        end
    end
end

% Create figure for plotting results
fig = figure('Position', [100, 100, 900, 600]);
plot(matrix_sizes, average_errors, '-o', 'LineWidth', 1.5);
xlabel('Matrix Size (n)', 'FontSize', 16);
ylabel('Average Relative Residual', 'FontSize', 16);
legend('Algo-LU','Algo-Frob','Algo-QR','Algo-SVD', 'Location', 'northwest', 'FontSize', 14);
title('Accuracy Comparison (Left Normwise Relative Residual)', 'FontSize', 18);
%title('Accuracy Comparison (Right Normwise Relative Residual)', 'FontSize', 18);
set(gca, 'FontSize', 14);
xlim([100 2000]); % Set x-axis limit for better visualization
ylim([0 3.5e-15]); % Set y-axis limit for better visualization
grid on;

% Display detected anomalies
fprintf('Detected anomalies:\n');
for i = 1:length(anomalies)
    a = anomalies{i};
    fprintf('Method %d | Size %dx%d | +%.2f%% over average | cond(A) = %.2e\n', ...
        a.method, a.matrix_size, a.matrix_size, a.percent_over, a.cond_num);
end

% Save the figure as high-resolution PNG
% filename = fullfile(pwd, 'Fig_5_2_2_part_2_1.png');
% exportgraphics(fig, filename, 'Resolution', 300);

% filename = fullfile(pwd, 'Fig_5_2_2_part_2_2.png');
% exportgraphics(fig, filename, 'Resolution', 300);