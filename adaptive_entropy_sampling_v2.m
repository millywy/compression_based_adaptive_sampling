% adaptive_entropy_sampling_v2.m
% ACC-driven three-level adaptive sampling atop WFPV.
% Uses entropy_controller.m to pick fs_next (LOW=6.25, MID=12.5, HIGH=25).

clear;

%% Configuration
fs0 = 125;              % original sampling rate
fs_hi = 25;
fs_mid = 12.5;
fs_lo = 12.5;
fs_proc = 25;           % internal WFPV processing rate
fs_acc = 25;            % fixed ACC control stream
FFTres = 1024;
WFlength = 15;
CutoffFreqHzBP = [0.4 4]; % bandpass for WFPV
CutoffFreqHzSearch = [1 3];
window_sec = 8;
step_sec = 2;
nbits_entropy = 2;
DEBUG_LOG = false;

IDData = {'DATA_01_TYPE01','DATA_02_TYPE02','DATA_03_TYPE02','DATA_04_TYPE02',...
    'DATA_05_TYPE02','DATA_06_TYPE02','DATA_07_TYPE02','DATA_08_TYPE02','DATA_09_TYPE02',...
    'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02','DATA_S04_T01',...
    'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
    'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
    'TEST_S07_T02', 'TEST_S08_T01'};

[b125, a125] = butter(4, CutoffFreqHzBP/(fs0/2), 'bandpass');

%% Metrics containers
logRec = struct([]);
myError = nan(1, numel(IDData));

%% Process each recording
for idnb = 1:numel(IDData)
    load(['Data/' IDData{idnb}], 'sig');
    if idnb > 13, ch = [1 2 3 4 5]; else, ch = [2 3 4 5 6]; end

    window = window_sec * fs0;
    step   = step_sec   * fs0;
    windowNb = floor((size(sig,2)-window)/step) + 1;
    if windowNb < 1, continue; end

    BPM_est = zeros(1, windowNb);
    FsUsed  = zeros(1, windowNb);
    Hacc    = zeros(1, windowNb);
    ACCmag25_log = zeros(1, windowNb);

    % Controller state (let ctrl_entropy_online initialize on first call)
    state_ctrl = [];
    fs_cur = fs_lo;

    % WFPV states
    state_hi = init_mode_state();
    state_mid = init_mode_state();
    state_lo = init_mode_state();

    % Pre-filter full recording at 125 Hz once for ACC entropy
    sig_sel = sig(ch, :);
    sig_filt = zeros(size(sig_sel));
    for c = 1:size(sig_sel,1)
        sig_filt(c,:) = filter(b125, a125, sig_sel(c,:));
    end

    for i = 1:windowNb
        curSegment = (i-1)*step+1 : (i-1)*step+window;
        curDataRaw = sig_sel(:, curSegment);      % raw for WFPV (125 Hz)
        curDataFilt = sig_filt(:, curSegment);    % filtered for ACC entropy

        % ACC control stream at fixed 25 Hz
        curAcc25 = do_resample(curDataFilt(3:5, :), fs0, fs_acc);
        ACCmag25 = sqrt(curAcc25(1,:).^2 + curAcc25(2,:).^2 + curAcc25(3,:).^2);
        ACCmag25_log(i) = mean(ACCmag25);
        Hacc(i) = entropy_proxy_arith_o1(ACCmag25, nbits_entropy);

        % Select WFPV state based on current fs_cur
        if abs(fs_cur - fs_hi) < eps
            state_in = state_hi;
        elseif abs(fs_cur - fs_mid) < eps
            state_in = state_mid;
        else
            state_in = state_lo;
        end

        % Run WFPV frame (handles resample to fs_cur then to fs_proc)
        [BPM_est(i), state_out] = wfpv_one_frame(curDataRaw, fs0, fs_cur, fs_proc, FFTres, WFlength, CutoffFreqHzSearch, state_in, i, BPM_est, idnb);

        % Save back mode state
        if abs(fs_cur - fs_hi) < eps
            state_hi = state_out;
        elseif abs(fs_cur - fs_mid) < eps
            state_mid = state_out;
        else
            state_lo = state_out;
        end

        FsUsed(i) = fs_cur;

        % Decide next fs using entropy_controller (3-level)
        accMag_norm = (ACCmag25_log(i)-min(ACCmag25_log(max(1,i-30):i)))/max(eps, max(ACCmag25_log(max(1,i-30):i))-min(ACCmag25_log(max(1,i-30):i)));
        [fs_next, state_ctrl] = ctrl_entropy_2state(Hacc(i), accMag_norm, 2, state_ctrl);
        fs_cur = fs_next;

        if DEBUG_LOG
            fprintf('rec %02d win %03d mode=%s fs=%4.1f Hacc=%.3f accMag=%.3f\n', idnb, i, state_ctrl.mode, fs_cur, Hacc(i), accMag_norm);
        end
    end

    % Ground truth and error
    if idnb > 13
        load(['Data/True' IDData{idnb}(5:end)], 'BPM0');
    else
        load(['Data/' IDData{idnb} '_BPMtrace'], 'BPM0');
    end
    frames = min(length(BPM_est), length(BPM0));
    myError(idnb) = mean(abs(BPM0(1:frames) - BPM_est(1:frames)'));

    % Store log
    logRec(idnb).ID = IDData{idnb};
    logRec(idnb).FsUsed = FsUsed;
    logRec(idnb).Hacc = Hacc;
    logRec(idnb).ACC_mag_25 = ACCmag25_log(1:frames);
    logRec(idnb).BPM_est = BPM_est(1:frames);
    logRec(idnb).BPM0 = BPM0(1:frames);
    logRec(idnb).entropyStats = struct( ...
        'median', median(Hacc), ...
        'p25', prctile(Hacc,25), ...
        'p75', prctile(Hacc,75), ...
        'max', max(Hacc));
end

%% Aggregate metrics
MAE_all = mean(myError, 'omitnan');
MAE_train = mean(myError(1:12), 'omitnan');
MAE_test = mean(myError(13:end), 'omitnan');
fprintf('\n=== Adaptive Entropy Sampling v2 Results ===\n');
fprintf('Err12=%2.2f, Err11=%2.2f, ErrAll=%2.2f\n', MAE_train, MAE_test, MAE_all);

%% Helpers
function state = init_mode_state()
state.W1_FFTi = [];
state.W11_FFTi = [];
state.W2_FFTi = [];
state.W21_FFTi = [];
state.prevFFT = [];
state.rangeIdx = [];
state.FreqRange = [];
end

function [BPM_val, state] = wfpv_one_frame(curDataRaw, fs0, fs_adc, fs_proc, FFTres, WFlength, searchHz, state, i, BPM_est, idnb)
% Bandpass at 125 Hz, resample to fs_adc, then to 25 Hz internal, then run original WFPV logic
[b,a] = butter(4, [0.4 4]/(fs0/2), 'bandpass');
curDataFilt = zeros(size(curDataRaw));
for c = 1:size(curDataRaw,1)
    curDataFilt(c,:) = filter(b,a,curDataRaw(c,:));
end
% resample to fs_adc
curData_adc = do_resample(curDataFilt, fs0, fs_adc);
% resample to internal 25 Hz if needed
if abs(fs_adc - fs_proc) < eps
    curData = curData_adc; fs = fs_proc;
else
    curData = do_resample(curData_adc, fs_adc, fs_proc);
    fs = fs_proc;
end

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

% Phase vocoder
FreqRangePPG = FreqRange;
if ~isempty(state.prevFFT) && length(state.prevFFT)==length(PPG_ave_FFT)
    for ii = 1:numel(FreqRangePPG)
        curPhase = angle(PPG_ave_FFT(ii));
        prevPhase = angle(state.prevFFT(ii));
        vocoder = zeros(1, 20);
        for n = 1:20
            vocoder(n) = ((curPhase - prevPhase) + (2 * pi * (n - 1))) / (2 * pi * 2);
        end
        [~, deltaidx] = min(abs(vocoder - FreqRange(ii)));
        FreqRangePPG(ii) = vocoder(deltaidx);
    end
end
FreqRangePPG = moving(FreqRangePPG, 3);
state.prevFFT = PPG_ave_FFT;

% Wiener filtering history within mode
WC1 = WFlength; WC2 = WFlength;

state.W1_FFTi = [state.W1_FFTi; (abs(PPG_ave_FFT))/max(abs(PPG_ave_FFT)+eps)];
W1_hist = state.W1_FFTi(max(1,end-WC1+1):end, :);
W1_PPG_ave_FFT_ALL = mean(W1_hist,1);
W1_PPG_ave_FFT_ALL_norm = W1_PPG_ave_FFT_ALL / max(W1_PPG_ave_FFT_ALL + eps);
W1_ACC_X_FFT_norm = (abs(ACC_X_FFT)) / max(abs(ACC_X_FFT) + eps);
W1_ACC_Y_FFT_norm = (abs(ACC_Y_FFT)) / max(abs(ACC_Y_FFT) + eps);
W1_ACC_Z_FFT_norm = (abs(ACC_Z_FFT)) / max(abs(ACC_Z_FFT) + eps);
WF1 = (1 - 1/3 * (W1_ACC_X_FFT_norm + W1_ACC_Y_FFT_norm + W1_ACC_Z_FFT_norm) ./ (W1_PPG_ave_FFT_ALL_norm + eps));
WF1(WF1 < 0) = -1;
W1_PPG_ave_FFT_Clean = abs(PPG_ave_FFT) .* WF1;

state.W11_FFTi = [state.W11_FFTi; abs(PPG_ave_FFT).^2];
W11_hist = state.W11_FFTi(max(1,end-WC1+1):end, :);
W11_PPG_ave_FFT_ALL = mean(W11_hist,1);
W11_PPG_ave_FFT_ALL_norm = W11_PPG_ave_FFT_ALL / max(W11_PPG_ave_FFT_ALL + eps);
W11_ACC_X_FFT_norm = (abs(ACC_X_FFT) .^ 2) / max(abs(ACC_X_FFT .^ 2) + eps);
W11_ACC_Y_FFT_norm = (abs(ACC_Y_FFT) .^ 2) / max(abs(ACC_Y_FFT .^ 2) + eps);
W11_ACC_Z_FFT_norm = (abs(ACC_Z_FFT) .^ 2) / max(abs(ACC_Z_FFT.^ 2) + eps);
WF11 = (1 - 1/3 * (W11_ACC_X_FFT_norm + W11_ACC_Y_FFT_norm + W11_ACC_Z_FFT_norm) ./ (W11_PPG_ave_FFT_ALL_norm + eps));
W11_PPG_ave_FFT_Clean = abs(PPG_ave_FFT) .* WF11;

state.W2_FFTi = [state.W2_FFTi; (abs(PPG_ave_FFT))/max(abs(PPG_ave_FFT)+eps)];
W2_hist = state.W2_FFTi(max(1,end-WC2+1):end, :);
W2_PPG_ave_FFT_ALL = mean(W2_hist,1);
W2_PPG_ave_FFT_ALL_norm = W2_PPG_ave_FFT_ALL / max(W2_PPG_ave_FFT_ALL + eps);
W2_ACC_X_FFT_norm = (abs(ACC_X_FFT)) / max(abs(ACC_X_FFT) + eps);
W2_ACC_Y_FFT_norm = (abs(ACC_Y_FFT)) / max(abs(ACC_Y_FFT) + eps);
W2_ACC_Z_FFT_norm = (abs(ACC_Z_FFT)) / max(abs(ACC_Z_FFT) + eps);
WF2 = W2_PPG_ave_FFT_ALL_norm ./ (((W2_ACC_X_FFT_norm + W2_ACC_Y_FFT_norm + W2_ACC_Z_FFT_norm) / 3) + W2_PPG_ave_FFT_ALL_norm + eps);
W2_PPG_ave_FFT_Clean = abs(PPG_ave_FFT) .* WF2;

state.W21_FFTi = [state.W21_FFTi; abs(PPG_ave_FFT).^2];
W21_hist = state.W21_FFTi(max(1,end-WC2+1):end, :);
W21_PPG_ave_FFT_ALL = mean(W21_hist,1);
W21_PPG_ave_FFT_ALL_norm = W21_PPG_ave_FFT_ALL / max(W21_PPG_ave_FFT_ALL + eps);
W21_ACC_X_FFT_norm = abs(ACC_X_FFT) .^ 2 / max(abs(ACC_X_FFT) .^ 2 + eps);
W21_ACC_Y_FFT_norm = abs(ACC_Y_FFT) .^ 2 / max(abs(ACC_Y_FFT) .^ 2 + eps);
W21_ACC_Z_FFT_norm = abs(ACC_Z_FFT) .^ 2 / max(abs(ACC_Z_FFT) .^ 2 + eps);
WF21 = W21_PPG_ave_FFT_ALL_norm ./ (((W21_ACC_X_FFT_norm + W21_ACC_Y_FFT_norm + W21_ACC_Z_FFT_norm) / 3) + W21_PPG_ave_FFT_ALL_norm + eps);
W21_PPG_ave_FFT_Clean = abs(PPG_ave_FFT) .* WF21;

W1_PPG_ave_FFT_Clean = W1_PPG_ave_FFT_Clean / std(W1_PPG_ave_FFT_Clean + eps);
W11_PPG_ave_FFT_Clean = W11_PPG_ave_FFT_Clean / std(W11_PPG_ave_FFT_Clean + eps);
W2_PPG_ave_FFT_Clean = W2_PPG_ave_FFT_Clean / std(W2_PPG_ave_FFT_Clean + eps);
W21_PPG_ave_FFT_Clean = W21_PPG_ave_FFT_Clean / std(W21_PPG_ave_FFT_Clean + eps);

PPG_ave_FFT_FIN = W1_PPG_ave_FFT_Clean + W2_PPG_ave_FFT_Clean;

% History-tracked peak picking
hist_int = 25;
if idnb > 12
    if i > 15, hist_int = max(abs(diff(BPM_est(1:i-1)))) + 5; end
else
    if i > 30, hist_int = max(abs(diff(BPM_est(1:i-1)))) + 5; end
end

if isempty(state.rangeIdx)
    [~, idx] = max(PPG_ave_FFT_FIN);
    BPM_val = FreqRangePPG(idx(1)) * 60;
    state.rangeIdx = idx(1) - round(hist_int / ((FreqRange(2) - FreqRange(1)) * 60)) : ...
                     idx(1) + round(hist_int / ((FreqRange(2) - FreqRange(1)) * 60));
else
    [~, idx] = max(PPG_ave_FFT_FIN(state.rangeIdx));
    BPM_val = FreqRangePPG(state.rangeIdx(idx(1))) * 60;
    state.rangeIdx = state.rangeIdx(idx(1)) - round(hist_int / ((FreqRange(2) - FreqRange(1)) * 60)) : ...
                     state.rangeIdx(idx(1)) + round(hist_int / ((FreqRange(2) - FreqRange(1)) * 60));
end
state.rangeIdx(state.rangeIdx < 1) = [];
state.rangeIdx(state.rangeIdx > length(FreqRange)) = [];

if i > 5 && abs(BPM_val - BPM_est(i - 1)) > 5
    ddd = polyfit(1:length(BPM_est(max(1, i - 5):i - 1)), BPM_est(max(1, i - 5):i - 1), 1);
    BPM_val = 0.8 * BPM_val + 0.2 * polyval(ddd, length(BPM_est(max(1, i - 5):i - 1)) + 1);
end

mul = 0.1;
if i > 1
    BPM_val = BPM_val + sum(sign(BPM_est(max(2, i - 6):i) - BPM_est(max(1, i - 7):i - 1)) * mul);
end
end

function sig_res = do_resample(sig, fs_in, fs_out)
    if abs(fs_out - fs_in) < eps
        sig_res = sig; return;
    end
    ratio = fs_out / fs_in;
    [nCh, nSamp] = size(sig);
    if fs_out < fs_in && abs(fs_in/fs_out - round(fs_in/fs_out)) < 1e-9
        decim = round(fs_in/fs_out);
        sig_res = sig(:, 1:decim:end);
    elseif fs_out > fs_in && abs(ratio - round(ratio)) < 1e-9
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
% Aggregate metrics and plots (similar to v1)
MAE_all = mean(myError, 'omitnan');
MAE_train = mean(myError(1:12), 'omitnan');
MAE_test = mean(myError(13:end), 'omitnan');
fprintf('\n=== Adaptive Entropy Sampling v2 Results ===\n');
fprintf('Err12=%2.2f, Err11=%2.2f, ErrAll=%2.2f\n', MAE_train, MAE_test, MAE_all);

% Bland-Altman and correlation
fullBPM0 = [];
fullBPM = [];
for rr = 1:numel(logRec)
    fullBPM0 = [fullBPM0, logRec(rr).BPM0(:)'];
    fullBPM  = [fullBPM,  logRec(rr).BPM_est(:)'];
end
fprintf('Generating Bland-Altman plot...\n');
[~, figBA] = BlandAltman(fullBPM0', fullBPM', {'Ground truth HR','Estimated HR'});
if exist('sgtitle','file') && ~isempty(figBA), figure(figBA); sgtitle('Bland-Altman (Adaptive Entropy Sampling v2)'); end
tmp = corrcoef(fullBPM0, fullBPM);
fprintf('Overall correlation coefficient: %.4f\n', tmp(1,2));

% Mode usage stats
all_fs_used = [logRec.FsUsed];
pct_hi = 100 * mean(all_fs_used==fs_hi);
pct_lo = 100 * mean(all_fs_used==fs_lo);
fprintf('Mode usage: high=%.1f%% low=%.1f%%\n', pct_hi, pct_lo);

save('adaptive_entropy_logs_v2.mat', 'logRec', 'myError', 'MAE_all', 'MAE_train', 'MAE_test');

%% Selected recordings for comparison
selRecs = [2 10 16 20];
for idx = 1:numel(selRecs)
    r = selRecs(idx);
    if r <= numel(logRec) && ~isempty(logRec(r).BPM0)
        win_count = numel(logRec(r).Hacc);
        frames = 1:win_count;

        % HR comparison
        figure;
        plot(logRec(r).BPM0,'ro'); hold on; plot(logRec(r).BPM_est,'o','Color','blue');
        title(sprintf('Recording %d (Adaptive Entropy Sampling v2)', r));
        xlabel('Time (frames)'); ylabel('HR (BPM)'); legend({'Ground truth','Estimates'});
        grid on;

        % Combined overlay: Hacc, FsUsed, ACC magnitude (dHacc_raw not logged in v2)
        figure;
        ax1 = subplot(3,1,1); plot(frames, logRec(r).Hacc, '-'); ylabel('Hacc (bits)'); title(sprintf('Recording %d - Entropy & Sampling', r));
        ax2 = subplot(3,1,2); stairs(frames, logRec(r).FsUsed, '-'); ylabel('FsUsed (Hz)'); grid on;
        ax3 = subplot(3,1,3); plot(frames, logRec(r).ACC_mag_25, '-'); ylabel('ACC mag'); xlabel('Time (frames)'); grid on;
        linkaxes([ax1 ax2 ax3],'x'); grid(ax1,'on');
    end
end
