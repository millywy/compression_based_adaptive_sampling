%% Bland-Altman Analysis for HR Estimation
% This script demonstrates how to use BlandAltman.m to compare 
% HR estimations with ground truth

% Modify your main script to collect all estimates and ground truth
% Then run this analysis

%% Basic Usage
% BlandAltman(data1, data2)
% - Creates correlation plot (left) and Bland-Altman plot (right)
% - data1: Ground truth BPM values (vector)
% - data2: Estimated BPM values (vector, same length as data1)

%% Example 1: Simple Comparison
% BlandAltman(ground_truth_bpm, estimated_bpm);

%% Example 2: With Labels and Units
% BlandAltman(ground_truth_bpm, estimated_bpm, ...
%     {'Ground Truth', 'Estimated BPM', 'BPM'});
%     % labels: {data1_name, data2_name, units}

%% Example 3: Full Configuration
% BlandAltman(ground_truth_bpm, estimated_bpm, ...
%     {'Ground Truth', 'Estimated BPM', 'BPM'}, ...    % label
%     'HR Estimation Validation', ...                    % title
%     {'Training', 'Test'}, ...                          % group names (if comparing multiple groups)
%     {'eq'; 'r2'; 'SSE'; 'n'}, ...                      % correlation info
%     {'RPC(%)'; 'CV'}, ...                              % Bland-Altman info
%     'auto');                                            % axis limits

% BlandAltman(ground_truth, estimated_values, labels, title, groupnames, correlation_info, BA_info)

% Example 4: Customized Info Display
% BlandAltman(BPM_ground_truth, BPM_estimated, ...
%     {'Ground Truth', 'Estimated BPM', 'BPM'}, ...
%     'HR Estimation Validation', ...
%     {}, ...
%     {'eq'; 'r2'; 'n'}, ...
%     {'RPC(%)'; 'CV'});
%% Correlation Info Options (displayed on correlation plot):
% 'eq'     - Slope and intercept equation (y = ax + b)
% 'r'      - Pearson correlation coefficient
% 'r2'     - Pearson R-squared
% 'rho'    - Spearman rank correlation
% 'SSE'    - Sum of squared errors
% 'n'      - Number of data points
% Default: {'eq'; 'r2'; 'SSE'; 'n'}

%% Bland-Altman Info Options (displayed on BA plot):
% 'RPC'    - Reproducibility coefficient (1.96*SD)
% 'RPC(%)' - RPC as percentage of mean
% 'CV'     - Coefficient of variation (%)
% 'p'      - P-value from t-test
% Default: {'RPC(%)'; 'CV'}

%% Return Values:
% cr = reproducibility coefficient (1.96*SD of differences)
% fig = figure handles (contains both subplots)
% sstruct = structure with statistics:
%   - N: number of data points
%   - CR: reproducibility coefficient
%   - r: Pearson correlation
%   - r2: R-squared
%   - SSE: sum of squared errors
%   - rho: Spearman correlation
%   - Slope: regression slope
%   - Intercept: regression intercept

%% How to Modify PPG_HRest_offline_ICASSP2017.m

% Step 1: Collect all BPM estimates and ground truth
% Add this BEFORE the end of the main loop (after computing myErrorN):

% global_BPM0 = [];       % Store all ground truth
% global_BPM_est_N = [];  % Store all estimates

% Then in the main loop, after computing BPM_est_N:
% global_BPM0 = [global_BPM0; BPM0(:)];
% global_BPM_est_N = [global_BPM_est_N; BPM_est_N(:)];

% Step 2: Create Bland-Altman plot after the main loop:
% [cr, fig, stats] = BlandAltman(global_BPM0, global_BPM_est_N, ...
%     {'Ground Truth', 'Viterbi Estimate', 'BPM'}, ...
%     'PPG Heart Rate Estimation - Bland-Altman Analysis', ...
%     {}, ...
%     {'eq'; 'r2'; 'rho'; 'n'}, ...
%     {'RPC(%)'; 'CV'}, ...
%     'auto');

% fprintf('\n=== Bland-Altman Statistics ===\n');
% fprintf('Reproducibility Coefficient: %.2f BPM\n', cr);
% fprintf('Pearson r: %.4f\n', stats.r);
% fprintf('R-squared: %.4f\n', stats.r2);
% fprintf('Slope: %.4f, Intercept: %.4f\n', stats.Slope, stats.Intercept);

%% Output Interpretation:

% LEFT PLOT (Correlation):
% - X-axis: Ground truth BPM
% - Y-axis: Estimated BPM
% - Black line: Linear regression fit
% - Grey dashed line: Perfect agreement (y = x)
% - Points closer to grey line = better agreement
% - r2 > 0.95: Excellent correlation
% - r2 > 0.80: Good correlation

% RIGHT PLOT (Bland-Altman):
% - X-axis: Mean of ground truth and estimate
% - Y-axis: Difference (Estimate - Ground Truth)
% - Black line: Mean difference
% - Dashed lines: ±1.96*SD (limits of agreement)
% - Points should be randomly scattered around zero
% - Narrower band (small ±1.96*SD) = more consistent estimates
% - RPC(%) < 10%: Excellent agreement
% - RPC(%) < 20%: Good agreement

%% Example Output Values:
% A good HR estimation algorithm should show:
% - r2 > 0.90
% - RPC < 10 BPM (or RPC < 15% for variable HR)
% - Points randomly scattered around zero on BA plot
% - No systematic bias (mean difference ≈ 0)


% Key Statistics Returned:
% RPC (Reproducibility Coefficient): 1.96×SD of differences
% r² (R-squared): Correlation strength (0-1)
% CV (Coefficient of Variation): Relative error (%)
% Slope/Intercept: Regression fit parameters
