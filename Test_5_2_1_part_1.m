% This program compares the computational speed of 
% Algo-LU, Algo-Frob, Algo-QR, Algo-SVD
% The comparison is performed on random Gaussian matrices to evaluate
% speed performance across different problem dimensions

% Get the current project directory and add helper functions to the path
projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'help_functions'));

% Define matrix sizes to test (from 100 to 2000 in steps of 100)
matrix_sizes = 100:100:2000;  % 20 different sizes in total
num_trials = 10;              % Number of trials for each matrix size
num_methods = 4;              % Number of methods to compare
average_times = zeros(length(matrix_sizes), num_methods); % Store average times

anomalies = {};  % Initialize cell array to record anomalies: method, matrix, time, average time, deviation, condition number

for idx = 1:length(matrix_sizes)
    n = matrix_sizes(idx);
    times = zeros(num_trials, num_methods); % Store times for each trial

    % Perform multiple trials for statistical reliability
    for trial = 1:num_trials
        rng(trial);  % Set fixed random seed for reproducibility
        M = randn(n) + 1i * randn(n);  % Generate random Gaussian complex matrix

        % Algo-LU
        tic;
        X1 = inverse_via_LU(M);
        times(trial, 1) = toc;

        % Algo-Frob
        tic;
        X2 = inverse_Frobenius(M);
        times(trial, 2) = toc;

        % Algo-QR
        tic;
        X3 = inverse_via_QR(M);
        times(trial, 3) = toc;

        % Algo-SVD
        tic;
        X4 = inverse_via_SVD(M);
        times(trial, 4) = toc;
    end

    % Calculate average times for each method
    avg_time = mean(times, 1);
    average_times(idx, :) = avg_time;

    % Check for anomalies (trials where time exceeds 150% of average)
    for method = 1:num_methods
        for trial = 1:num_trials
            t = times(trial, method);
            if t > 1.5 * avg_time(method)
                condM = cond(M);
                anomalies{end+1} = struct( ...
                    'matrix_size', n, ...
                    'method', method, ...
                    'trial_time', t, ...
                    'avg_time', avg_time(method), ...
                    'percent_over', 100 * (t - avg_time(method)) / avg_time(method), ...
                    'cond_num', condM ...
                );
            end
        end
    end
end

% Create figure for plotting results
fig = figure('Position', [100, 100, 1100, 600]);
plot(matrix_sizes, average_times, '-o', 'LineWidth', 1.5);
xlabel('Matrix Size (n)', 'FontSize', 16);
ylabel('Average Runtime (s)', 'FontSize', 16);
legend('Algo-LU','Algo-Frob','Algo-QR','Algo-SVD', 'Location', 'northwest', 'FontSize', 14);
title('Speed Comparison (Average Runtime)', 'FontSize', 18);
set(gca, 'FontSize', 14);
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
% filename = fullfile(pwd, 'Fig_5_2_1_part_1.png');
% exportgraphics(fig, filename, 'Resolution', 300);