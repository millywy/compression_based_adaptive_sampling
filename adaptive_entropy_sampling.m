% adaptive_entropy_sampling.m
% Entropy-driven adaptive sampling atop the WFPV (TBME2017) pipeline.
% Two modes: fs_hi=25 Hz, fs_lo=12.5 Hz. Raw data assumed at fs0=125 Hz.

clear;  % close all;

%% Configuration
fs0 = 125;
fs_hi = 25;
fs_lo = 12.5;
fs_acc = 25;                 % fixed-rate control stream for ACC
FFTres = 1024;
WFlength = 15;                % Wiener averaging length (frames)
CutoffFreqHzBP = [0.4 4];     % bandpass at 125 Hz before decimation
CutoffFreqHzSearch = [1 3];   % HR search band (Hz)
window_sec = 8;
step_sec = 2;

% Adaptive entropy thresholds / hysteresis (tunable)
W_hist = 30;                  % history length in windows (~60s at 2s hop)
W_min = 10;                   % warm-up windows before normal adaptation
nbits_entropy = 4;            % quantization for entropy proxy
hi_hold_init = 3;             % min windows to stay high after switching up
ema_beta = 0.85;              % smoothing for entropy (ACC-based control)
k_up_lo = 3.0;                % sensitivity when in LOW mode
k_up_hi = 3.0;                % stricter when in HIGH mode
k_dn    = 1.0;                % stability requirement for downshift
N_up_level = 2;               % consecutive level-high windows to go HIGH
N_up_jump  = 1;               % consecutive jump triggers to go HIGH (burst on sudden change)
N_down  = 4;                  % consecutive stable windows to go LOW
cooldown_init = 3;            % block re-up for this many windows after going LOW
DEBUG_LOG = false;            % set true to print per-window debug

IDData = {'DATA_01_TYPE01','DATA_02_TYPE02','DATA_03_TYPE02','DATA_04_TYPE02',...
    'DATA_05_TYPE02','DATA_06_TYPE02','DATA_07_TYPE02','DATA_08_TYPE02','DATA_09_TYPE02',...
    'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02','DATA_S04_T01',...
    'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
    'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
    'TEST_S07_T02', 'TEST_S08_T01'};

%% Filter at original rate (reuse per recording)
[b125, a125] = butter(4, CutoffFreqHzBP/(fs0/2), 'bandpass');

%% Metrics containers
results = struct('fs', [], 'MAE_all', [], 'MAE_train', [], 'MAE_test', []);
logRec = struct([]);
myError = nan(1, numel(IDData));

%% Process each recording
for idnb = 1:numel(IDData)
    % Load data
    load(['Data/' IDData{idnb}], 'sig');
    if idnb > 13
        ch = [1 2 3 4 5];
    else
        ch = [2 3 4 5 6];
    end

    window = window_sec * fs0;
    step   = step_sec   * fs0;
    windowNb = floor((size(sig,2)-window)/step) + 1;

    BPM_est = zeros(1, windowNb);
    FsUsed  = zeros(1, windowNb);
    Hacc    = zeros(1, windowNb);      % ACC entropy at fixed 25 Hz
    Hacc_s  = zeros(1, windowNb);      % smoothed ACC entropy
    dHacc_raw = zeros(1, windowNb);    % raw ACC entropy diff (jump metric)
    Th_high_log = zeros(1, windowNb);  % store enter-high gate
    Th_low_log  = zeros(1, windowNb);  % store exit-low gate
    dH_up_log   = zeros(1, windowNb);  % jump threshold
    dH_dn_log   = zeros(1, windowNb);  % stability threshold

    % Mode states (keep WF history per mode)
    state_hi = init_mode_state();
    state_lo = init_mode_state();
    fs_cur = fs_hi;   % start high
    hi_hold = 0;
    state_mode = "LOW";
    cooldown = 0; up_count_level = 0; up_count_jump = 0; down_count = 0;
    switches_up = 0; switches_down = 0;

    for i = 1:windowNb
        curSegment = (i-1)*step+1 : (i-1)*step+window;
        curDataRaw = sig(ch, curSegment);

        % filter at 125 Hz
        curDataFilt = zeros(size(curDataRaw));
        for c = 1:size(curDataRaw,1)
            curDataFilt(c,:) = filter(b125, a125, curDataRaw(c,:));
        end

        % ACC control stream at fixed 25 Hz
        acc_decim = round(fs0/fs_acc); % should be 5
        curAcc25 = curDataFilt(3:5, 1:acc_decim:end);
        ACCmag25 = sqrt(curAcc25(1,:).^2 + curAcc25(2,:).^2 + curAcc25(3,:).^2);

        % PPG/ACC for WFPV at current mode
        if fs_cur == fs_hi
            decim = round(fs0/fs_hi); % 5
            curData = curDataFilt(:, 1:decim:end); % downsample along time (cols) to 25 Hz
            state = state_hi;
        else
            decim = round(fs0/fs_lo); % 10
            curData = curDataFilt(:, 1:decim:end); % downsample along time (cols) to 12.5 Hz
            state = state_lo;
        end

        % Run one-frame WFPV at fs_cur
        [BPM_est(i), state] = wfpv_one_frame(curData, fs_cur, FFTres, WFlength, CutoffFreqHzSearch, state, i, BPM_est, idnb);

        % Save back mode state
        if fs_cur == fs_hi
            state_hi = state;
        else
            state_lo = state;
        end

        % Entropy metrics: ACC control stream at fixed 25 Hz
        Hacc(i) = entropy_proxy_context(ACCmag25, nbits_entropy);
        if i==1
            Hacc_s(i) = Hacc(i);
            dHacc_raw(i) = 0;
        else
            Hacc_s(i) = ema_beta * Hacc_s(i-1) + (1-ema_beta) * Hacc(i);
            dHacc_raw(i) = Hacc(i) - Hacc(i-1); % jump metric (raw)
        end

        % Adaptive thresholds from recent history (ACC-based)
        if i < W_min
            Th_enter = inf; Th_exit = -inf; dH_dn = inf; dH_up = inf; sigma_d = 0;
            hi_hold = hi_hold_init; % warm-up: stay high
            up_count_level = 0; up_count_jump = 0; down_count = 0; cooldown = 0;
        else
            i0 = max(1, i-W_hist+1);
            H_hist = Hacc_s(i0:i);
            medH = median(H_hist);
            iqrH = prctile(H_hist,75) - prctile(H_hist,25) + eps;
            Th_enter = medH + 0.8 * iqrH; % gate to go HIGH
            Th_exit  = medH + 0.2 * iqrH; % gate to go LOW (hysteresis)

            dH_hist = abs(dHacc_raw(max(2,i-W_hist+1):i));
            if isempty(dH_hist), dH_hist = 0; end
            trim_thr = prctile(dH_hist,90);
            d_trim = dH_hist(dH_hist <= trim_thr);
            if isempty(d_trim), d_trim = 0; end
            sigma_d = 1.4826 * mad(d_trim,1) + 1e-4;
            k_up = k_up_lo;
            if state_mode == "HIGH"
                k_up = k_up_hi;
            end
            dH_up = k_up * sigma_d;
            dH_dn = k_dn  * sigma_d;
        end

        Th_high_log(i) = Th_enter; % store enter-high gate here
        Th_low_log(i)  = Th_exit;
        dH_up_log(i)   = dH_up;
        dH_dn_log(i)   = dH_dn;

        trig_up_jump = false; trig_up_level = false; stable = false; big_dip = false;

        % State machine controller (burst-mode)
        if i < W_min
            state_mode = "LOW";
            fs_cur = fs_lo;
        else
            switch state_mode
                case "LOW"
                    fs_cur = fs_lo;
                    if cooldown > 0
                        cooldown = cooldown - 1;
                        up_count_level = 0; up_count_jump = 0;
                    else
                        trig_up_jump = dHacc_raw(i) > dH_up;
                        trig_up_level = Hacc_s(i) > Th_enter;
                        if trig_up_level
                            up_count_level = up_count_level + 1;
                        else
                            up_count_level = 0;
                        end
                        if trig_up_jump
                            up_count_jump = up_count_jump + 1;
                        else
                            up_count_jump = 0;
                        end
                        if (up_count_level >= N_up_level) || (up_count_jump >= N_up_jump)
                            state_mode = "HIGH";
                            fs_cur = fs_hi;
                            hi_hold = hi_hold_init;
                            up_count_level = 0; up_count_jump = 0; down_count = 0;
                            switches_up = switches_up + 1;
                        end
                    end
                case "HIGH"
                    fs_cur = fs_hi;
                    big_dip = dHacc_raw(i) < -dH_up; % immediate down on large negative jump
                    if hi_hold > 0
                        hi_hold = hi_hold - 1;
                    else
                        stable = (Hacc_s(i) < Th_exit) && (abs(dHacc_raw(i)) < k_dn*sigma_d);
                        if stable
                            down_count = down_count + 1;
                        else
                            down_count = 0;
                        end
                        if down_count >= N_down
                            state_mode = "LOW";
                            fs_cur = fs_lo;
                            cooldown = cooldown_init;
                            up_count_level = 0; up_count_jump = 0; down_count = 0;
                            switches_down = switches_down + 1;
                        end
                    end
            end
        end
        FsUsed(i) = fs_cur;

        % Optional debug (toggle with DEBUG_LOG)
        if DEBUG_LOG
            fprintf(['rec %02d win %03d mode=%s fs=%4.1f Hacc_s=%.3f dHraw=%.3f ', ...
                'ThEnter=%.3f ThExit=%.3f sig=%.3f jump=%d level=%d dip=%d stable=%d upL=%d upJ=%d dn=%d cool=%d hold=%d\n'], ...
                idnb, i, state_mode, fs_cur, Hacc_s(i), dHacc_raw(i), Th_enter, Th_exit, sigma_d, ...
                trig_up_jump, trig_up_level, big_dip, stable, up_count_level, up_count_jump, down_count, cooldown, hi_hold);
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
    logRec(idnb).Hacc_s = Hacc_s;
    logRec(idnb).dHacc_raw = dHacc_raw;
    logRec(idnb).BPM_est = BPM_est(1:frames);
    logRec(idnb).BPM0 = BPM0(1:frames);
    logRec(idnb).Th_high_log = Th_high_log;
    logRec(idnb).Th_low_log  = Th_low_log;
    logRec(idnb).dH_up_log   = dH_up_log;
    logRec(idnb).dH_dn_log   = dH_dn_log;
    logRec(idnb).switches_up = switches_up;
    logRec(idnb).switches_down = switches_down;

    % Entropy summary stats
    stats.Hacc.median = median(Hacc);
    stats.Hacc.p25 = prctile(Hacc,25);
    stats.Hacc.p75 = prctile(Hacc,75);
    stats.Hacc.max = max(Hacc);
    logRec(idnb).entropyStats = stats;
end

%% Aggregate metrics
MAE_all = mean(myError, 'omitnan');
MAE_train = mean(myError(1:12), 'omitnan');
MAE_test = mean(myError(13:end), 'omitnan');
fprintf('\n=== Adaptive Entropy Sampling Results ===\n');
fprintf('Err12=%2.2f, Err11=%2.2f, ErrAll=%2.2f\n', MAE_train, MAE_test, MAE_all);
fprintf('Individual recording errors (BPM):\n');
fprintf(' ');
fprintf('%4.2f ', myError);
fprintf('\n');

% Bland-Altman and correlation
fullBPM0 = [];
fullBPM = [];
for rr = 1:numel(logRec)
    fullBPM0 = [fullBPM0, logRec(rr).BPM0(:)'];
    fullBPM  = [fullBPM,  logRec(rr).BPM_est(:)'];
end
fprintf('Generating Bland-Altman plot...\n');
[~, figBA] = BlandAltman(fullBPM0', fullBPM', {'Ground truth HR','Estimated HR'});
if exist('sgtitle','file') && ~isempty(figBA), figure(figBA); sgtitle('Bland-Altman (Adaptive Entropy Sampling)'); end
tmp = corrcoef(fullBPM0, fullBPM);
fprintf('Overall correlation coefficient: %.4f\n', tmp(1,2));

% Mode usage stats
all_fs_used = [logRec.FsUsed];
pct_hi = 100 * mean(all_fs_used==fs_hi);
pct_lo = 100 * mean(all_fs_used==fs_lo);
fprintf('Mode usage: high=%.1f%% low=%.1f%%\n', pct_hi, pct_lo);

save('adaptive_entropy_logs.mat', 'logRec', 'myError', 'MAE_all', 'MAE_train', 'MAE_test');

%% Selected recordings for comparison
selRecs = [6 12 20];
for idx = 1:numel(selRecs)
    r = selRecs(idx);
    if r <= numel(logRec) && ~isempty(logRec(r).BPM0)
        win_count = numel(logRec(r).Hacc);
        t_sec = (0:win_count-1) * step_sec; % window start times
        frames = 1:win_count;

        % HR comparison
        figure;
        plot(logRec(r).BPM0,'ro'); hold on; plot(logRec(r).BPM_est,'o','Color','blue');
        title(sprintf('Recording %d (Adaptive Entropy Sampling)', r));
        xlabel('Time (frames)'); ylabel('HR (BPM)'); legend({'Ground truth','Estimates'});
        grid on;

        % Hacc and FsUsed over time (time axis in seconds, optional frame ticks on top)
        figure;
        yyaxis left; plot(frames, logRec(r).Hacc, '-'); ylabel('Hacc (bits)');
        yyaxis right; stairs(frames, logRec(r).FsUsed, '-'); ylabel('FsUsed (Hz)');
        xlabel('Time (frames)'); title(sprintf('Entropy & FsUsed - Recording %d', r));
        grid on;
    end
end

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

function [BPM_val, state] = wfpv_one_frame(curData, fs, FFTres, WFlength, searchHz, state, i, BPM_est, idnb)
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
W11_ACC_Z_FFT_norm = (abs(ACC_Z_FFT) .^ 2) / max(abs(ACC_Z_FFT .^ 2) + eps);
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

% History tracking
hist_int = 25;
if idnb > 12
    if i > 15, hist_int = max(abs(diff(BPM_est(max(1,i-15):i-1)))) + 5; end
else
    if i > 30, hist_int = max(abs(diff(BPM_est(max(1,i-30):i-1)))) + 5; end
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
    prev_seq = BPM_est(max(1,i-7):i-1);
    if numel(prev_seq) > 1
        BPM_val = BPM_val + sum(sign(diff(prev_seq))) * mul;
    end
end
end


%% compressors choices for entropy 

%% Entropy proxy via delta symbols and histogram
function H = entropy_proxy(x, nbits)
x = (x - mean(x)) / (std(x) + eps);
x = max(min(x, 3), -3);
levels = 2^nbits;
q = floor((x - min(x)) / (max(x) - min(x) + eps) * (levels - 1));
if numel(q) < 2
    H = 0;
    return;
end
dq = diff(q);
[counts, ~] = histcounts(dq, -levels:levels);
p = counts / sum(counts);
p = p(p>0);
H = -sum(p .* log2(p));
end

%% ENTROPY_PROXY_CONTEXT Entropy estimator via context modeling (order-1) + ideal arithmetic length.
function H = entropy_proxy_context(x, nbits)
%ENTROPY_PROXY_CONTEXT Entropy estimator via context modeling (order-1) + ideal arithmetic length.
%   H = entropy_proxy_context(x, nbits)
%   - x: 1D signal
%   - nbits: number of quantization bits
%   Returns: estimated bits per delta-symbol (lower = more predictable/compressible)

    x = x(:)';
    if numel(x) < 4
        H = 0;
        return;
    end

    % 1) Normalize + clip (robust-ish)
    x = (x - mean(x)) / (std(x) + eps);
    x = max(min(x, 3), -3);

    % 2) Uniform quantization to L levels
    L = 2^nbits;
    % map [-3,3] -> [0, L-1]
    q = round((x + 3) * (L-1) / 6);
    q = min(max(q, 0), L-1);

    % 3) Delta symbols (centered), then shift to nonnegative alphabet
    d = diff(q);                       % range approx [-(L-1), +(L-1)]
    S = 2*(L-1) + 1;                   % alphabet size
    sym = d + (L-1) + 1;               % -> [1..S]

    % 4) Order-1 context model: P(sym_t | sym_{t-1})
    % We'll compute the ideal arithmetic codelength: sum -log2 P
    alpha = 1;                         % Laplace smoothing
    counts = zeros(S, S);              % counts(prev, curr)

    % prime with first symbol as "prev"
    prev = sym(1);

    bits = 0;
    for t = 2:numel(sym)
        curr = sym(t);

        row = counts(prev, :);
        denom = sum(row) + alpha*S;
        numer = row(curr) + alpha;
        p = numer / denom;

        bits = bits - log2(p);

        % update model
        counts(prev, curr) = counts(prev, curr) + 1;
        prev = curr;
    end

    H = bits / (numel(sym)-1);         % bits per symbol (delta)
end
