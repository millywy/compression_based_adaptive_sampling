% baseline_sampling_sweep.m
% Sweep fixed sampling rates for the WFPV (Wiener Filter + Phase Vocoder) TBME 2017 pipeline.
% Uses the existing logic from PPG_WFPV_TBME2017, but runs the full pipeline
% at multiple uniform sampling frequencies without modifying the original script.

clear;  % close all;

%% Configuration
fs0 = 125;  % original sampling rate of the dataset
fs_list = [125 100 50 25 20 15 12.5 10];
FFTres = 1024;
FFTres_base = 1024;  % tuned at 25 Hz
fs_base = 25;
WFlength = 15;          % Wiener averaging length (frames)
CutoffFreqHzBP = [0.4 4];  % bandpass for PPG/ACC filtering (Hz)
CutoffFreqHzSearch = [1 3]; % search band for HR peaks (Hz, 60–180 BPM)

IDData = {'DATA_01_TYPE01','DATA_02_TYPE02','DATA_03_TYPE02','DATA_04_TYPE02',...
    'DATA_05_TYPE02','DATA_06_TYPE02','DATA_07_TYPE02','DATA_08_TYPE02','DATA_09_TYPE02',...
    'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02','DATA_S04_T01',...
    'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
    'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
    'TEST_S07_T02', 'TEST_S08_T01'};

%% Sweep over sampling rates
results = struct('fs', [], 'MAE_all', [], 'MAE_train', [], 'MAE_test', []);
for k = 1:numel(fs_list)
    fs = fs_list(k);
    fprintf('\n=== fs_target = %.2f Hz ===\n', fs);
    myError = nan(1, numel(IDData));

    for idnb = 1:numel(IDData)
        % Load data
        load(['Data/' IDData{idnb}], 'sig');
        if idnb > 13
            ch = [1 2 3 4 5];
        else
            ch = [2 3 4 5 6];
        end

        % Resample to target rate with anti-aliasing
        sig_res = resample_to_fs(sig(ch, :), fs0, fs);

        % Run WFPV pipeline at this fs
        BPM_est = run_wfpv_record(sig_res, fs, fs_base, FFTres_base, WFlength, CutoffFreqHzBP, CutoffFreqHzSearch, idnb);

        % Load ground truth BPM trace
        if idnb > 13
            load(['Data/True' IDData{idnb}(5:end)], 'BPM0');
        else
            load(['Data/' IDData{idnb} '_BPMtrace'], 'BPM0');
        end

        frames = min(length(BPM_est), length(BPM0));
        myError(idnb) = mean(abs(BPM0(1:frames) - BPM_est(1:frames)'));
    end

    results(k).fs = fs;
    results(k).MAE_all = mean(myError, 'omitnan');
    results(k).MAE_train = mean(myError(1:12), 'omitnan');
    results(k).MAE_test = mean(myError(13:end), 'omitnan');

    fprintf('MAE all: %.2f | train: %.2f | test: %.2f\n', ...
        results(k).MAE_all, results(k).MAE_train, results(k).MAE_test);
end

%% Plot MAE vs sampling frequency
figure;
plot([results.fs], [results.MAE_all], '-o', 'LineWidth', 1.25);
xlabel('Sampling frequency (Hz)');
ylabel('MAE (BPM)');
title('WFPV MAE vs sampling frequency');
grid on;

%% Helper: resample with anti-aliasing
function sig_res = resample_to_fs(sig, fs_in, fs_out)
    % sig: channels x samples
    if abs(fs_out - fs_in) < eps
        sig_res = sig;
        return;
    end

    r = fs_out / fs_in;
    [nCh, nSamp] = size(sig);
    expected_len = round(nSamp * r);

    % Integer decimation when possible
    if fs_out < fs_in && abs(fs_in / fs_out - round(fs_in / fs_out)) < 1e-6
        decim = round(fs_in / fs_out);
        sig_res = zeros(nCh, ceil(nSamp / decim));
        for c = 1:nCh
            sig_res(c, :) = decimate(sig(c, :), decim);
        end
    else
        sig_res = zeros(nCh, expected_len);
        for c = 1:nCh
            resampled = resample(sig(c, :), fs_out, fs_in);
            % Trim or pad to expected length
            if length(resampled) > expected_len
                sig_res(c, :) = resampled(1:expected_len);
            else
                sig_res(c, 1:length(resampled)) = resampled;
            end
        end
    end
end

%% Helper: WFPV pipeline for one recording at a given fs
function BPM_est = run_wfpv_record(sig, fs, fs_base, FFTres_base, WFlength, bpHz, searchHz, idnb)
    FFTres = 2^nextpow2(round(FFTres_base * (fs/fs_base)));
    [b, a] = butter(4, bpHz / (fs / 2), 'bandpass');
    window = round(8 * fs);
    step = round(2 * fs);

    windowNb = floor((size(sig, 2) - window) / step) + 1;
    if windowNb < 1
        BPM_est = [];
        return;
    end

    BPM_est = zeros(1, windowNb);
    rangeIdx = [];
    clear W1_FFTi W11_FFTi W2_FFTi W21_FFTi W1_PPG_ave_FFT_Clean W2_PPG_ave_FFT_Clean W11_PPG_ave_FFT_Clean PPG_ave_FFT_FIN W21_PPG_ave_FFT_Clean PPG_ave_FFT_FIN;

    for i = 1:windowNb
        curSegment = (i - 1) * step + 1 : (i - 1) * step + window;
        curData = sig(:, curSegment);

        PPG1 = filter(b, a, curData(1, :));
        PPG2 = filter(b, a, curData(2, :));
        ACC_X = filter(b, a, curData(3, :));
        ACC_Y = filter(b, a, curData(4, :));
        ACC_Z = filter(b, a, curData(5, :));
        PPG_ave = 0.5 * (PPG1 - mean(PPG1)) / (std(PPG1) + eps) + 0.5 * (PPG2 - mean(PPG2)) / (std(PPG2) + eps);

        % Periodogram
        PPG_ave_FFT = fft(PPG_ave, FFTres);
        FreqRange = linspace(0, fs, size(PPG_ave_FFT, 2));

        [~, lowR] = min(abs(FreqRange - searchHz(1)));
        [~, highR] = min(abs(FreqRange - searchHz(2)));
        FreqRange = FreqRange(lowR:highR);
        PPG_ave_FFT = PPG_ave_FFT(lowR:highR);
        ACC_X_FFT = fft(ACC_X, FFTres); ACC_X_FFT = ACC_X_FFT(lowR:highR);
        ACC_Y_FFT = fft(ACC_Y, FFTres); ACC_Y_FFT = ACC_Y_FFT(lowR:highR);
        ACC_Z_FFT = fft(ACC_Z, FFTres); ACC_Z_FFT = ACC_Z_FFT(lowR:highR);

        % Phase vocoder refinement
        FreqRangePPG = FreqRange;
        if i > 1
            for ii = 1:numel(FreqRangePPG)
                curPhase = angle(PPG_ave_FFT(ii));
                prevPhase = angle(PPG_ave_FFTpr(ii));
                vocoder = zeros(1, 20);
                for n = 1:20
                    vocoder(n) = ((curPhase - prevPhase) + (2 * pi * (n - 1))) / (2 * pi * 2);
                end
                [~, deltaidx] = min(abs(vocoder - FreqRange(ii)));
                FreqRangePPG(ii) = vocoder(deltaidx);
            end
        end
        FreqRangePPG = moving(FreqRangePPG, 3);
        PPG_ave_FFTpr = PPG_ave_FFT;

        % Wiener filtering
        WC1 = WFlength; WC2 = WFlength;

        W1_FFTi(i, :) = (abs(PPG_ave_FFT)) / max(abs(PPG_ave_FFT) + eps);
        if i == 1, W1_PPG_ave_FFT_ALL = W1_FFTi(i, :); else, W1_PPG_ave_FFT_ALL = mean(W1_FFTi(max(1, i - WC1):i, :), 1); end
        W1_PPG_ave_FFT_ALL_norm = W1_PPG_ave_FFT_ALL / max(W1_PPG_ave_FFT_ALL + eps);
        W1_ACC_X_FFT_norm = (abs(ACC_X_FFT)) / max(abs(ACC_X_FFT) + eps);
        W1_ACC_Y_FFT_norm = (abs(ACC_Y_FFT)) / max(abs(ACC_Y_FFT) + eps);
        W1_ACC_Z_FFT_norm = (abs(ACC_Z_FFT)) / max(abs(ACC_Z_FFT) + eps);
        WF1 = (1 - 1/3 * (W1_ACC_X_FFT_norm + W1_ACC_Y_FFT_norm + W1_ACC_Z_FFT_norm) ./ (W1_PPG_ave_FFT_ALL_norm + eps));
        WF1(WF1 < 0) = -1;
        W1_PPG_ave_FFT_Clean(i, :) = abs(PPG_ave_FFT) .* WF1;

        W11_FFTi(i, :) = abs(PPG_ave_FFT) .^ 2;
        if i == 1, W11_PPG_ave_FFT_ALL = W11_FFTi(i, :); else, W11_PPG_ave_FFT_ALL = mean(W11_FFTi(max(1, i - WC1):i, :), 1); end
        W11_PPG_ave_FFT_ALL_norm = W11_PPG_ave_FFT_ALL / max(W11_PPG_ave_FFT_ALL + eps);
        W11_ACC_X_FFT_norm = (abs(ACC_X_FFT) .^ 2) / max(abs(ACC_X_FFT .^ 2) + eps);
        W11_ACC_Y_FFT_norm = (abs(ACC_Y_FFT) .^ 2) / max(abs(ACC_Y_FFT .^ 2) + eps);
        W11_ACC_Z_FFT_norm = (abs(ACC_Z_FFT) .^ 2) / max(abs(ACC_Z_FFT .^ 2) + eps);
        WF11 = (1 - 1/3 * (W11_ACC_X_FFT_norm + W11_ACC_Y_FFT_norm + W11_ACC_Z_FFT_norm) ./ (W11_PPG_ave_FFT_ALL_norm + eps));
        W11_PPG_ave_FFT_Clean(i, :) = abs(PPG_ave_FFT) .* WF11;

        W2_FFTi(i, :) = (abs(PPG_ave_FFT)) / max(abs(PPG_ave_FFT) + eps);
        if i == 1, W2_PPG_ave_FFT_ALL = W2_FFTi(i, :); else, W2_PPG_ave_FFT_ALL = mean(W2_FFTi(max(1, i - WC2):i, :), 1); end
        W2_PPG_ave_FFT_ALL_norm = W2_PPG_ave_FFT_ALL / max(W2_PPG_ave_FFT_ALL + eps);
        W2_ACC_X_FFT_norm = (abs(ACC_X_FFT)) / max(abs(ACC_X_FFT) + eps);
        W2_ACC_Y_FFT_norm = (abs(ACC_Y_FFT)) / max(abs(ACC_Y_FFT) + eps);
        W2_ACC_Z_FFT_norm = (abs(ACC_Z_FFT)) / max(abs(ACC_Z_FFT) + eps);
        WF2 = W2_PPG_ave_FFT_ALL_norm ./ (((W2_ACC_X_FFT_norm + W2_ACC_Y_FFT_norm + W2_ACC_Z_FFT_norm) / 3) + W2_PPG_ave_FFT_ALL_norm + eps);
        W2_PPG_ave_FFT_Clean(i, :) = abs(PPG_ave_FFT) .* WF2;
        W2_FFTi(i, :) = (W2_PPG_ave_FFT_Clean(i, :)) / max(W2_PPG_ave_FFT_Clean(i, :) + eps);

        W21_FFTi(i, :) = abs(PPG_ave_FFT) .^ 2;
        if i == 1, W21_PPG_ave_FFT_ALL = W21_FFTi(i, :); else, W21_PPG_ave_FFT_ALL = mean(W21_FFTi(max(1, i - WC2):i, :), 1); end
        W21_PPG_ave_FFT_ALL_norm = W21_PPG_ave_FFT_ALL / max(W21_PPG_ave_FFT_ALL + eps);
        W21_ACC_X_FFT_norm = abs(ACC_X_FFT) .^ 2 / max(abs(ACC_X_FFT) .^ 2 + eps);
        W21_ACC_Y_FFT_norm = abs(ACC_Y_FFT) .^ 2 / max(abs(ACC_Y_FFT) .^ 2 + eps);
        W21_ACC_Z_FFT_norm = abs(ACC_Z_FFT) .^ 2 / max(abs(ACC_Z_FFT) .^ 2 + eps);
        WF21 = W21_PPG_ave_FFT_ALL_norm ./ (((W21_ACC_X_FFT_norm + W21_ACC_Y_FFT_norm + W21_ACC_Z_FFT_norm) / 3) + W21_PPG_ave_FFT_ALL_norm + eps);
        W21_PPG_ave_FFT_Clean(i, :) = abs(PPG_ave_FFT) .* WF21;
        W21_FFTi(i, :) = W21_PPG_ave_FFT_Clean(i, :) .^ 2;

        W1_PPG_ave_FFT_Clean(i, :) = W1_PPG_ave_FFT_Clean(i, :) / std(W1_PPG_ave_FFT_Clean(i, :) + eps);
        W11_PPG_ave_FFT_Clean(i, :) = W11_PPG_ave_FFT_Clean(i, :) / std(W11_PPG_ave_FFT_Clean(i, :) + eps);
        W2_PPG_ave_FFT_Clean(i, :) = W2_PPG_ave_FFT_Clean(i, :) / std(W2_PPG_ave_FFT_Clean(i, :) + eps);
        W21_PPG_ave_FFT_Clean(i, :) = W21_PPG_ave_FFT_Clean(i, :) / std(W21_PPG_ave_FFT_Clean(i, :) + eps);

        PPG_ave_FFT_FIN(i, :) = W1_PPG_ave_FFT_Clean(i, :) + W2_PPG_ave_FFT_Clean(i, :);

        % History-tracked peak picking
        hist_int = 25;
        if idnb > 12
            if i > 15, hist_int = max(abs(diff(BPM_est(1:i-1)))) + 5; end
        else
            if i > 30, hist_int = max(abs(diff(BPM_est(1:i-1)))) + 5; end
        end

        if isempty(rangeIdx)
            [~, idx] = max(PPG_ave_FFT_FIN(i, :));
            BPM_est(i) = FreqRangePPG(idx(1)) * 60;
            rangeIdx = idx(1) - round(hist_int / ((FreqRange(2) - FreqRange(1)) * 60)) : ...
                       idx(1) + round(hist_int / ((FreqRange(2) - FreqRange(1)) * 60));
        else
            [~, idx] = max(PPG_ave_FFT_FIN(i, rangeIdx));
            BPM_est(i) = FreqRangePPG(rangeIdx(idx(1))) * 60;
            rangeIdx = rangeIdx(idx(1)) - round(hist_int / ((FreqRange(2) - FreqRange(1)) * 60)) : ...
                       rangeIdx(idx(1)) + round(hist_int / ((FreqRange(2) - FreqRange(1)) * 60));
        end
        rangeIdx(rangeIdx < 1) = [];
        rangeIdx(rangeIdx > length(FreqRange)) = [];

        if i > 5 && abs(BPM_est(i) - BPM_est(i - 1)) > 5
            ddd = polyfit(1:length(BPM_est(max(1, i - 5):i - 1)), BPM_est(max(1, i - 5):i - 1), 1);
            BPM_est(i) = 0.8 * BPM_est(i) + 0.2 * polyval(ddd, length(BPM_est(max(1, i - 5):i - 1)) + 1);
        end

        mul = 0.1;
        BPM_est(i) = BPM_est(i) + sum(sign(BPM_est(max(2, i - 6):i) - BPM_est(max(1, i - 7):i - 1)) * mul);
    end
end
