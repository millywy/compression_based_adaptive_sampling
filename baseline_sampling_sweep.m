% baseline_sampling_sweep.m
% Sweep ADC sampling rates while keeping the internal WFPV pipeline fixed at 25 Hz
% (matching PPG_WFPV_TBME2017: filter at 125 Hz -> window -> resample to ADC -> resample to 25 Hz -> WFPV).

clear;  % close all;

%% Configuration
fs0 = 125;  % original sampling rate of the dataset
fs_list = [25 12.5 6.25];
%fs_list = [25 20 12.5 10 6.25 5];
fs_proc = 25;           % fixed internal rate for WFPV
FFTres = 1024;          % FFT length at 25 Hz (original)
WFlength = 15;          % Wiener averaging length (frames)
CutoffFreqHzBP = [0.4 4];   % bandpass for PPG/ACC filtering (Hz) at 125 Hz
CutoffFreqHzSearch = [1 3]; % search band for HR peaks (Hz, 60–180 BPM) at 25 Hz

IDData = {'DATA_01_TYPE01','DATA_02_TYPE02','DATA_03_TYPE02','DATA_04_TYPE02',...
    'DATA_05_TYPE02','DATA_06_TYPE02','DATA_07_TYPE02','DATA_08_TYPE02','DATA_09_TYPE02',...
    'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02','DATA_S04_T01',...
    'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
    'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
    'TEST_S07_T02', 'TEST_S08_T01'};
traceLog = repmat(struct('ID',[],'BPM0',[],'est',struct('fs',[],'BPM',[])),1,numel(IDData));

%% Sweep over sampling rates
results = struct('fs', [], 'MAE_all', [], 'MAE_train', [], 'MAE_test', []);
for k = 1:numel(fs_list)
    fs_adc = fs_list(k);
    fprintf('\n=== fs_target (ADC) = %.2f Hz ===\n', fs_adc);
    myError = nan(1, numel(IDData));
    fullBPM0 = [];
    fullBPM = [];

    for idnb = 1:numel(IDData)
        % Load data
        load(['Data/' IDData{idnb}], 'sig');
        if idnb > 13
            ch = [1 2 3 4 5];
        else
            ch = [2 3 4 5 6];
        end

        % Run WFPV pipeline with per-window resampling to fs_adc then to fs_proc (clean helper)
        % Run WFPV pipeline with entropy helper (returns Hacc for plotting)
        [BPM_est, Hacc_est] = run_wfpv_entropy(sig(ch, :), fs0, fs_adc, fs_proc, FFTres, WFlength, CutoffFreqHzSearch, idnb, CutoffFreqHzBP);

        % Load ground truth BPM trace
        if idnb > 13
            load(['Data/True' IDData{idnb}(5:end)], 'BPM0');
        else
            load(['Data/' IDData{idnb} '_BPMtrace'], 'BPM0');
        end

        frames = min(length(BPM_est), length(BPM0));
        myError(idnb) = mean(abs(BPM0(1:frames) - BPM_est(1:frames)'));
        fullBPM0 = [fullBPM0, BPM0(1:frames)'];
        fullBPM  = [fullBPM, BPM_est(1:frames)];

        % Cache traces for multi-fs plotting after sweep
        if isempty(traceLog(idnb).ID)
            traceLog(idnb).ID = IDData{idnb};
            traceLog(idnb).BPM0 = BPM0(:)';
            traceLog(idnb).est = struct('fs', {}, 'BPM', {});
        end
        estEntry = struct('fs', fs_adc, 'BPM', BPM_est(1:frames), 'Hacc', Hacc_est(1:frames));
        if isempty(traceLog(idnb).est)
            traceLog(idnb).est = estEntry;
        else
            traceLog(idnb).est(end+1) = estEntry;
        end

        % Plot selected recordings 9 and 14 for comparison
        % if idnb==2 || idnb==10 || idnb==16 || idnb==20
        %     figure;
        %     plot(BPM0,'ro'); hold on; plot(BPM_est(1:frames),'o','Color','blue');
        %     title(sprintf('Recording %d at fs=%.2f Hz', idnb, fs_adc));
        %     xlabel('Time (frames)'); ylabel('HR (BPM)'); legend({'Ground truth','Estimates'});
        % end
    end

    results(k).fs = fs_adc;
    results(k).MAE_all = mean(myError, 'omitnan');
    results(k).MAE_train = mean(myError(1:12), 'omitnan');
    results(k).MAE_test = mean(myError(13:end), 'omitnan');

    fprintf('MAE all: %.2f | train: %.2f | test: %.2f\n', ...
        results(k).MAE_all, results(k).MAE_train, results(k).MAE_test);

    fprintf('Generating Bland-Altman plot for fs=%.2f...\n', fs_adc);
    [~, figBA] = BlandAltman(fullBPM0', fullBPM', {'Ground truth HR','Estimated HR'});
    if exist('sgtitle','file') && ~isempty(figBA), figure(figBA); sgtitle(sprintf('Bland-Altman (fs=%.2f Hz)', fs_adc)); end
    tmp = corrcoef(fullBPM0, fullBPM);
    fprintf('Overall correlation coefficient (fs=%.2f): %.4f\n', fs_adc, tmp(1,2));
end

%% Plot selected recordings with all fs overlays
selRecs = [1 2 3 4 5 6 7 8 9 10];
for idx = 1:numel(selRecs)
    r = selRecs(idx);
    if r > numel(traceLog) || isempty(traceLog(r).ID) || isempty(traceLog(r).est)
        continue;
    end
    % sort estimates by fs for consistent legend
    est_fs = [traceLog(r).est.fs];
    [~, ord] = sort(est_fs);
    est_sorted = traceLog(r).est(ord);

    figure;
    % plot ground truth (trim to min length across estimates)
    min_len = min(cellfun(@(x) numel(x), [{traceLog(r).BPM0}, arrayfun(@(e) e.BPM, est_sorted, 'UniformOutput', false)]));
    hold on;
    plot(traceLog(r).BPM0(1:min_len), 'o', 'Color','black', 'LineWidth', 1.25, 'DisplayName', 'Ground truth');
    cmap = lines(numel(est_sorted));
    for ee = 1:numel(est_sorted)
        frames = min(min_len, numel(est_sorted(ee).BPM));
        plot(est_sorted(ee).BPM(1:frames), 'o', 'Color', cmap(ee,:), 'LineWidth', 1.0, ...
            'DisplayName', sprintf('Estimate fs=%.2f', est_sorted(ee).fs));
    end
    title(sprintf('Recording %d: BPM estimates vs ground truth', r));
    xlabel('Time (frames)'); ylabel('HR (BPM)');
    legend('Location','best'); grid on;

    % Hacc overlay on same subplot (right axis) for quick comparison
    yyaxis right;
    min_len_hacc = min(cellfun(@(x) numel(x), arrayfun(@(e) e.Hacc, est_sorted, 'UniformOutput', false)));
    for ee = 1:numel(est_sorted)
        frames = min(min_len_hacc, numel(est_sorted(ee).Hacc));
        plot(est_sorted(ee).Hacc(1:frames), '--', 'Color', cmap(ee,:), 'LineWidth', 1.0, ...
            'DisplayName', sprintf('Hacc fs=%.2f', est_sorted(ee).fs));
    end
    ylabel('Hacc (bits)');
    yyaxis left;
end

%% Plot MAE vs sampling frequency
figure;
plot([results.fs], [results.MAE_all], '-o', 'LineWidth', 1.25);
xlabel('Sampling frequency (Hz)');
ylabel('MAE (BPM)');
title('WFPV MAE vs sampling frequency (25 Hz internal)');
grid on;
