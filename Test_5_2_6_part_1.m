% This program compares the computational speed of 
% Algo-LU, Algo-Frob, Algo-Chol, Algo-Frob-Chol
% The comparison is performed on random HPD matrices to evaluate
% speed performance across different problem dimensions

projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'help_functions'));

matrix_sizes = 100:100:2000;
num_trials = 10;
num_methods = 4;
average_times = zeros(length(matrix_sizes), num_methods);

anomalies = {};

for idx = 1:length(matrix_sizes)
    n = matrix_sizes(idx);
    times = zeros(num_trials, num_methods);

    for trial = 1:num_trials
        rng(trial);
        N = randn(n) + 1i * randn(n);
        M = N' * N + 0.01 * eye(n);

        % Algo-LU
        tic;
        X1 = inverse_via_LU(M);
        times(trial, 1) = toc;

        % Algo-Frob
        tic;
        X2 = inverse_Frobenius(M);
        times(trial, 2) = toc;

        % Algo-Chol
        tic;
        X3 = inverse_via_Chol(M);
        times(trial, 3) = toc;

        % Algo-Frob-Chol
        tic;
        X4 = inverse_Frobenius_Chol(M);
        times(trial, 4) = toc;
    end

    avg_time = mean(times, 1);
    average_times(idx, :) = avg_time;

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

fig = figure('Position', [100, 100, 1100, 600]);
plot(matrix_sizes, average_times, '-o', 'LineWidth', 1.5);
xlabel('Matrix Size (n)', 'FontSize', 16);
ylabel('Average Runtime (s)', 'FontSize', 16);
legend('Algo-LU','Algo-Frob','Algo-Chol','Algo-Frob-Chol', 'Location', 'northwest', 'FontSize', 14);
title('Speed Comparison (Average Runtime)', 'FontSize', 18);
set(gca, 'FontSize', 14);
xlim([100 2000]);
grid on;

fprintf('Detected anomalies:\n');
for i = 1:length(anomalies)
    a = anomalies{i};
    fprintf('Method %d | Size %dx%d | +%.2f%% over average | cond(A) = %.2e\n', ...
        a.method, a.matrix_size, a.matrix_size, a.percent_over, a.cond_num);
end

% filename = fullfile(pwd, 'Fig_5_2_1_part_1.png');
% exportgraphics(fig, filename, 'Resolution', 300);