%% Helper: WFPV pipeline for one recording at fixed 25 Hz
function [BPM_est, Hacc, ACCmag25] = run_wfpv_record(sig_raw, fs0, fs_adc, fs_proc, FFTres, WFlength, searchHz, idnb, bpHz)
    % Bandpass at 125 Hz, then per-window resample to fs_adc, then (if needed) back to 25 Hz for the original WFPV logic
    [b125,a125] = butter(4, bpHz/(fs0/2), 'bandpass');
    fs_acc = 25;           % fixed ACC control stream for entropy
    nbits_entropy = 2;     % quantization for entropy proxy
    window125 = round(8 * fs0);
    step125 = round(2 * fs0);

    windowNb = floor((size(sig_raw, 2) - window125) / step125) + 1;
    if windowNb < 1
        BPM_est = [];
        Hacc = [];
        ACCmag25= [];
        return;
    end

    BPM_est = zeros(1, windowNb);
    Hacc = zeros(1, windowNb);
    % dHacc = zeros(1, windowNb);
    ACCmag25= zeros(1, windowNb);
    

    rangeIdx = [];
    clear W1_FFTi W11_FFTi W2_FFTi W21_FFTi W1_PPG_ave_FFT_Clean W2_PPG_ave_FFT_Clean W11_PPG_ave_FFT_Clean PPG_ave_FFT_FIN W21_PPG_ave_FFT_Clean PPG_ave_FFT_FIN;


    for i = 1:windowNb
        curSegment = (i - 1) * step125 + 1 : (i - 1) * step125 + window125;
        curDataRaw = sig_raw(:, curSegment);

        % filter at 125 Hz
        curDataFilt = zeros(size(curDataRaw));
        for c = 1:size(curDataRaw,1)
            curDataFilt(c,:) = filter(b125,a125,curDataRaw(c,:));
        end

        % resample to fs_adc (ADC rate) using decimate when integer ratio
        % resample to fs_adc (ADC rate) using decimate when integer ratio; else use rational resample
        curData_adc = do_resample(curDataFilt, fs0, fs_adc);
        if abs(fs_adc - fs_proc) < eps
            curData = curData_adc; fs = fs_proc;
        else
            curData = do_resample(curData_adc, fs_adc, fs_proc);
            fs = fs_proc;
        end

        % ACC control stream at fixed 25 Hz
        curAcc25 = do_resample(curDataFilt(3:5, :), fs0, fs_acc);
        accMagVec = sqrt(curAcc25(1,:).^2 + curAcc25(2,:).^2 + curAcc25(3,:).^2);
        ACCmag25(i) = mean(accMagVec);
        Hacc(i) = entropy_proxy_arith_o1(accMagVec, nbits_entropy);

        PPG1 = curData(1, :);
        PPG2 = curData(2, :);
        ACC_X = curData(3, :);
        ACC_Y = curData(4, :);
        ACC_Z = curData(5, :);
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

function sig_res = do_resample(sig, fs_in, fs_out)
    % Resample helper that uses decimation when integer ratio, otherwise rational resample.
    if abs(fs_out - fs_in) < eps
        sig_res = sig;
        return;
    end
    ratio = fs_out / fs_in;
    [nCh, nSamp] = size(sig);
    if fs_out < fs_in && abs(fs_in/fs_out - round(fs_in/fs_out)) < 1e-9
        decim = round(fs_in/fs_out);
        sig_res = sig(:, 1:decim:end);
    elseif fs_out > fs_in && abs(ratio - round(ratio)) < 1e-9
        % simple upsample by integer factor (zero-order hold via resample)
        up = round(ratio);
        expected_len = nSamp * up;
        sig_res = zeros(nCh, expected_len);
        for c = 1:nCh
            sig_res(c,:) = resample(sig(c,:), up, 1);
        end
    else
        [p,q] = rat(ratio,1e-6);
        expected_len = round(nSamp * p / q);
        sig_res = zeros(nCh, expected_len);
        for c = 1:nCh
            tmp = resample(sig(c,:), p, q);
            if length(tmp) > expected_len
                sig_res(c,:) = tmp(1:expected_len);
            else
                sig_res(c,1:length(tmp)) = tmp;
            end
        end
    end
end
