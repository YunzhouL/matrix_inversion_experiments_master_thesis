% This program is an alternative of Test_5_2_6_part_2_2
% The comparison is performed on random HPD matrices to evaluate
% relative forward error across different problem dimensions

projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'help_functions'));

matrix_sizes = 100:100:1000;
num_trials = 10;
num_methods = 4;
average_errors = zeros(length(matrix_sizes), num_methods);

anomalies = {};

for idx = 1:length(matrix_sizes)
    n = matrix_sizes(idx);
    errors = zeros(num_trials, num_methods);

    for trial = 1:num_trials
        rng(trial); 

        % Generate random HPD matrices
        N = randn(n) + 1i * randn(n);
        M = N' * N + 0.01 * eye(n);
        norm_M = norm(M, 'fro');

        X1 = inverse_via_LU(M);
        errors(trial, 1) = norm(X1 * M - eye(n), 'fro') / (norm(X1, 'fro') * norm_M);

        X2 = inverse_Frobenius(M);
        errors(trial, 2) = norm(X2 * M - eye(n), 'fro') / (norm(X2, 'fro') * norm_M);

        X3 = inverse_via_Chol(M);
        errors(trial, 3) = norm(X3 * M - eye(n), 'fro') / (norm(X3, 'fro') * norm_M);

        X4 = inverse_Frobenius_Chol(M);
        errors(trial, 4) = norm(X4 * M - eye(n), 'fro') / (norm(X4, 'fro') * norm_M);
    end

    avg_error = mean(errors, 1);
    average_errors(idx, :) = avg_error;

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

fig = figure('Position', [100, 100, 1100, 600]);
plot(matrix_sizes, average_errors, '-o', 'LineWidth', 1.5);
xlabel('Matrix Size (n)', 'FontSize', 16);
ylabel('Average Relative Forward Error (log scale)', 'FontSize', 16);
legend('Algo-LU','Algo-Frob','Algo-Chol','Algo-Frob-Chol','Location', 'northwest', 'FontSize', 14);
title('Accuracy Comparison ( M = N*N+0.01I )', 'FontSize', 18);
set(gca, 'FontSize', 14);
set(gca, 'YScale', 'log');
xlim([100 2000]);
grid on;

fprintf('Detected anomalies:\n');
for i = 1:length(anomalies)
    a = anomalies{i};
    fprintf('Method %d | Size %dx%d | +%.2f%% over average | cond(A) = %.2e\n', ...
        a.method, a.matrix_size, a.matrix_size, a.percent_over, a.cond_num);
end

% filename = fullfile(pwd, 'Fig_5_2_6_part_2.png');
% exportgraphics(fig, filename, 'Resolution', 300);