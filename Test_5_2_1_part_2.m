% This program compares the accuracy of Algo-LU, Algo-Frob, Algo-QR, Algo-SVD
% The comparison is performed on random Gaussian matrices to evaluate
% left relative residual across different problem dimensions

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
    errors = zeros(num_trials, num_methods);% Store errors for each trial

    % Perform multiple trials for statistical reliability
    for trial = 1:num_trials
        rng(trial);  % Set fixed random seed for reproducibility
        M = randn(n) + 1i * randn(n); % Generate random Gaussian complex matrix

        norm_M = norm(M, 'fro');

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
fig = figure('Position', [100, 100, 1100, 600]);
plot(matrix_sizes, average_errors, '-o', 'LineWidth', 1.5);
xlabel('Matrix Size (n)', 'FontSize', 16);
ylabel('Average Relative Residual (log scale)', 'FontSize', 16);
legend('Algo-LU','Algo-Frob','Algo-QR','Algo-SVD','Location', 'northwest', 'FontSize', 14);
title('Accuracy Comparison (Left Normwise Relative Residual)', 'FontSize', 18);
set(gca, 'FontSize', 14);
set(gca, 'YScale', 'log');
xlim([100 2000]); % Set x-axis limit for better visualization
grid on;

% Display detected anomalies
fprintf('Detected anomalies:\n');
for i = 1:length(anomalies)
    a = anomalies{i};
    fprintf('Method %d | Size %dx%d | +%.2f%% over average | cond(A) = %.2e\n', ...
        a.method, a.matrix_size, a.matrix_size, a.percent_over, a.cond_num);
end

% Save the figure as high-resolution PNG
% filename = fullfile(pwd, 'Fig_5_2_1_part_2.png');
% exportgraphics(fig, filename, 'Resolution', 300);