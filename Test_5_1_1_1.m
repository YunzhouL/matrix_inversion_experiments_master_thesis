% This program reproduces the experiment from the reference paper:
% "Complex matrix inversion via real matrix inversions"
% GitHub repository: https://github.com/zhen06/Complex-MatrixInversion
% This implementation uses MATLAB's backslash operator (\) for matrix inversion

% Get the current project directory and add helper functions to the path
projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'help_functions'));

% Define matrix sizes to test (from 2500 to 3500 in steps of 100)
matrix_sizes = 2500:100:3500;  % 11 different sizes in total
num_trials = 10;               % Number of trials for each matrix size
num_methods = 2;               % Number of methods to compare
average_times = zeros(length(matrix_sizes), num_methods); % Store average times

anomalies = {}; % Initialize cell array to record anomalies: method, matrix, time, average time, deviation, condition number

for idx = 1:length(matrix_sizes)
    n = matrix_sizes(idx);
    times = zeros(num_trials, num_methods); % Store times for each trial

    % Perform multiple trials for statistical reliability
    for trial = 1:num_trials
        rng(trial);  % Set fixed random seed for reproducibility

        % Generate random complex matrix M = A + i*B
        A = rand(n,n);
        B = rand(n,n);
        M = A + 1j*B;

        % Method 1: Direct inversion using MATLAB's backslash operator
        tic;
        M_inv1 = M \ eye(n);
        times(trial, 1) = toc;

        % Method 2: Frobenius-Schur inversion using backslash operator
        tic;
        X_2 = A\B;
        X_3 = (A + B*X_2) \ eye(n);
        X_4 = X_2*X_3;
        Minv2 = X_3 - 1j*X_4;
        times(trial, 2) = toc;

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
fig = figure('Position', [100, 100, 900, 600]);
plot(matrix_sizes, average_times, '-o', 'LineWidth', 1.5);
xlabel('Matrix Size (n)', 'FontSize', 16);
ylabel('Average Runtime (s)', 'FontSize', 16);
legend('Algo-LU','Algo-Frob', 'Location', 'northwest', 'FontSize', 14);
title('Speed Comparison using inv(A) inversion', 'FontSize', 18);
set(gca, 'FontSize', 14);
ylim([0 2]); % Set y-axis limit for better visualization
grid on;

% Display detected anomalies
fprintf('Detected anomalies:\n');
for i = 1:length(anomalies)
    a = anomalies{i};
    fprintf('Method %d | Size %dx%d | +%.2f%% over average | cond(A) = %.2e\n', ...
        a.method, a.matrix_size, a.matrix_size, a.percent_over, a.cond_num);
end

% Save the figure as high-resolution PNG
% filename = fullfile(pwd, 'Fig_5_1_1_1.png');
% exportgraphics(fig, filename, 'Resolution', 300);