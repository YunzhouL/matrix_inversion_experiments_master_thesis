% This program compares the accuracy of Algo-LU, Algo-Frob
% The comparison is performed on specially constructed random matrices
% M = A + iB, where κ(M) is signifin=cantly larger than κ(A), κ(B), κ(A + BA⁻¹B)
% and evaluate the relative forward error across different problem dimensions

projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'help_functions'));

matrix_sizes = 100:100:1000;
num_trials = 10;
num_methods = 2;
average_errors = zeros(length(matrix_sizes), num_methods);

anomalies = {};

for idx = 1:length(matrix_sizes)
    n = matrix_sizes(idx);
    errors = zeros(num_trials, num_methods);

    for trial = 1:num_trials
        rng(trial);

        % Generate M = A +iB with κ(A) = κ(B) = 10
        % κ(A + BA⁻¹B) <= 10, κ(M) > 10000
        t = 0.9999;        % Parameter t>0
        A = gallery('randsvd', n, 1e1, 2);
        [Q, R] = qr(A);
        
        % Construct N
        N = zeros(n);
        for j = 1:n/2
            N(2*j-1, 2*j) = -t;
            N(2*j, 2*j-1) =  t;
        end
        
        % Construct B = Q N R
        B = Q * N * R;
        M = A + 1i * B;

        M_inv = inverse_via_SVD(M);
        norm_M_inv = norm(M_inv,'fro');

        X1 = inverse_via_LU(M);
        errors(trial, 1) = norm(X1 - M_inv, 'fro') / norm_M_inv;

        X2 = inverse_Frobenius(M);
        errors(trial, 2) = norm(X2 - M_inv, 'fro') / norm_M_inv;
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
ylabel('Average Relative Forward Error', 'FontSize', 16);
legend('Algo-LU','Algo-Frob','Location', 'northwest', 'FontSize', 14);
title('Accuracy Comparison ( k(M)>10000 )', 'FontSize', 18);
set(gca, 'FontSize', 14);
set(gca, 'YScale', 'log');
grid on;

fprintf('Detected anomalies:\n');
for i = 1:length(anomalies)
    a = anomalies{i};
    fprintf('Method %d | Size %dx%d | +%.2f%% over average | cond(A) = %.2e\n', ...
        a.method, a.matrix_size, a.matrix_size, a.percent_over, a.cond_num);
end

% filename = fullfile(pwd, 'Fig_5_2_5_part_2.png');
% exportgraphics(fig, filename, 'Resolution', 300);