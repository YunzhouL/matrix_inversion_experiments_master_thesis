% This program compares the accuracy of Algo-LU and Algo-Frob
% The comparison is performed on specially constructed random matrices M = A + iB where:
% κ(M), κ(B), κ(A + BA⁻¹B) are held constant while  κ(A) is varied
% This experimental setup isolates the effect of κ(A) variations on
% algorithm accuracy, enabling analysis of how this specific condition number
% influences the relative forward error of each inversion method

projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'help_functions'));

matrix_sizes = 100:100:2000;
num_trials = 10;
num_methods = 2;
average_errors = zeros(length(matrix_sizes), num_methods);

anomalies = {};

for idx = 1:length(matrix_sizes)
    n = matrix_sizes(idx);
    errors = zeros(num_trials, num_methods);

    for trial = 1:num_trials
        rng(trial);

        % Generate M = A +iB with κ(M) = 14, κ(B) = 100, κ(A + BA⁻¹B) = 20,
        % κ(A) = 100 or 10
        U1 = RandOrthMat(n);
        V1 = RandOrthMat(n);
        % κ(A) = 100
        % [S1, S2] = generate_big_A(n, 100);
        % κ(A) = 10
        [S1, S2] = generate_small_A(n, 100);
        A = U1 * S1 * V1';
        B = U1 * S2 * V1';
        s = 1 ./ diag(S1 + 1i * S2);
        S_inv = diag(s);

        M = A + 1i * B;
        M_inv = V1 * S_inv * U1';
        norm_M_inv = norm(M_inv, 'fro');

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

condM = cond(S1 + 1i * S2);
condReal = cond(S1);
condImag = cond(S2);
condComp = cond(S1 + S2 * S1_inv * S2);
fprintf('cond(M) = %.2e | cond(A) = %.2e | cond(B) = %.2e | cond(A+BA^{-1}B) = %.2e\n', ...
condM, condReal, condImag, condComp);

fig = figure('Position', [100, 100, 900, 600]);
plot(matrix_sizes, average_errors, '-o', 'LineWidth', 1.5);
xlabel('Matrix Size (n)', 'FontSize', 16);
ylabel('Average Relative Forward Error', 'FontSize', 16);
legend('Algo-LU','Algo-Frob', 'Location', 'northwest', 'FontSize', 14);
% title('Accuracy Comparison ( k(A) = 100 )', 'FontSize', 18);
title('Accuracy Comparison ( k(A) = 10 )', 'FontSize', 18);
set(gca, 'FontSize', 14);
xlim([100 2000]);
grid on;

fprintf('Detected anomalies:\n');
for i = 1:length(anomalies)
    a = anomalies{i};
    fprintf('Method %d | Size %dx%d | +%.2f%% over average | cond(M) = %.2e | cond(A) = %.2e | cond(B) = %.2e\n  ', ...
        a.method, a.matrix_size, a.matrix_size, a.percent_over, a.cond_num, a.cond_real, a.cond_imag);
end

% filename = fullfile(pwd, 'Fig_5_2_4_part_1_1.png');
% exportgraphics(fig, filename, 'Resolution', 300);

% filename = fullfile(pwd, 'Fig_5_2_4_part_1_2.png');
% exportgraphics(fig, filename, 'Resolution', 300);