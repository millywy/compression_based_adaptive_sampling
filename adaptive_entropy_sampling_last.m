% adaptive_entropy_sampling.m
% Entropy-driven adaptive sampling atop the WFPV (TBME2017) pipeline.
% Two modes: fs_hi=25 Hz, fs_lo=12.5 Hz. Raw data assumed at fs0=125 Hz.

clear;  % close all;

%% Configuration
fs0 = 125;
fs_hi = 25;
fs_lo = 6.25;
fs_proc = 25;                 % internal WFPV processing rate (matches original)
FORCE_HIGH = false;            % set true to sanity-check fixed 25 Hz mode
fs_acc = 25;                 % fixed-rate control stream for ACC
FFTres = 1024;
WFlength = 15;          % Wiener averaging length (frames)
CutoffFreqHzBP = [0.4 4];     % bandpass at 125 Hz before decimation
CutoffFreqHzSearch = [1 3];   % HR search band (Hz)
window_sec = 8;
step_sec = 2;
fs_next = fs_lo;   % initially low

% Adaptive entropy thresholds / hysteresis (tunable)
% W_min = 10;                   % warm-up windows before normal adaptation
nbits_entropy = 2;            % quantization for entropy proxy
hi_hold = 15;                 % min windows to stay high after going HIGH
Th_hi = 0.15;                  % unstable = above this
N_look_back = 7;               % windows to look back for stable entropy
N_unstable = 3;                % unstable windows to go HIGH
Th_low = 0.13;                 % stable = below this
N_stable  = 5;                % min stable windows in high before going LOW

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

    sig_raw = sig(ch, :);

    window = window_sec * fs0;
    step   = step_sec   * fs0;
    windowNb = floor((size(sig_raw,2)-window)/step) + 1;

    if windowNb < 1
        BPM_est = [];
        Hacc = [];
        ACCmag25= [];
        return;
    end

    BPM_est = zeros(1, windowNb);
    FsUsed  = zeros(1, windowNb);
    Hacc    = zeros(1, windowNb);      % ACC entropy at fixed 25 Hz
    dHacc = zeros(1, windowNb);         % ACC entropy diff (jump metric)
    ACCmag25_log = zeros(1, windowNb); % mean ACC magnitude per window

    % Mode states (keep WF history per mode)
    state_hi = init_mode_state();
    state_lo = init_mode_state();
    state_in = init_mode_state();
    state_out = init_mode_state();
    fs_cur = fs_lo;   % start low
    state_mode = "LOW";  % start low

    rangeIdx = [];
    clear W1_FFTi W11_FFTi W2_FFTi W21_FFTi W1_PPG_ave_FFT_Clean W2_PPG_ave_FFT_Clean W11_PPG_ave_FFT_Clean PPG_ave_FFT_FIN W21_PPG_ave_FFT_Clean PPG_ave_FFT_FIN;

    
    for i = 1:windowNb
        curSegment = (i-1)*step+1 : (i-1)*step+window;
        curDataRaw = sig_raw(:, curSegment);

        % Run one-frame WFPV with internal resampling (match baseline helper)
        % load in correct mode state
        if fs_cur == fs_hi
            state_in = state_hi;
        else
            state_in = state_lo;
        end

        [BPM_est(i), state_out] = wfpv_one_frame_last(curDataRaw, fs0, fs_cur, fs_proc, FFTres, WFlength, CutoffFreqHzSearch, state_in, i, BPM_est, idnb, CutoffFreqHzBP);

        % Save back mode state
        if fs_cur == fs_hi
            state_hi = state_out;
        else
            state_lo = state_out;
        end

        % filter at 125 Hz
        curDataFilt2 = zeros(size(curDataRaw));
        for c = 1:size(curDataRaw,1)
            curDataFilt2(c,:) = filter(b125, a125, curDataRaw(c,:));
        end

        % ACC control stream at fixed 25 Hz
        curAcc_resampled = do_resample_last(curDataFilt2(3:5, :), fs0, fs_acc);
        ACCmag25 = sqrt(curAcc_resampled(1,:).^2 + curAcc_resampled(2,:).^2 + curAcc_resampled(3,:).^2);
        ACCmag25_log(i) = mean(ACCmag25);
        Hacc(i) = entropy_proxy_hist(ACCmag25, nbits_entropy);
        if i == 1
            dHacc(i) = 0;
        else
            dHacc(i) = Hacc(i) - Hacc(i-1);
        end

        % Adaptive sampling controller logic
        i0 = max(1, i - N_look_back + 1);
        recent = dHacc(i0:i);
        enough_hist = (i >= N_look_back);

        trig_up = false; trig_down = false; unstable_at_lo = false; stable_at_hi = false; 

        % State machine controller
        if FORCE_HIGH
            state_mode = "HIGH"; fs_next = fs_hi;
        else
            switch state_mode
                case "LOW"
                    fs_cur = fs_lo;
                    if enough_hist
                        n_bad = sum(abs(recent) > Th_hi);         % count unstable in last 5
                        if n_bad >= N_unstable
                            state_mode = "HIGH";
                            hi_timer = hi_hold;              % force HIGH for at least 15 windows
                            fs_next = fs_hi;
                        end
                    end
                case "HIGH"
                    fs_cur = fs_hi;
                    % Decrement hold timer while in HIGH
                    if hi_timer > 0
                        hi_timer = hi_timer - 1;
                    end

                    % After hold expires, exit only if ALL last 5 are stable (< Th_low)
                    if enough_hist && hi_timer == 0
                        n_good = sum(abs(recent) < Th_low);   % count stable in last N_look_back
                        if n_good >= N_stable
                            state_mode = "LOW";
                            fs_next = fs_lo;
                        end
                    end
                
            end
        end
        FsUsed(i) = fs_next;
        fs_cur = fs_next;

        % % Optional debug (toggle with DEBUG_LOG)
        % if DEBUG_LOG
        %     fprintf(['rec %02d win %03d mode=%s fs=%4.1f Hacc_s=%.3f dHraw=%.3f ', ...
        %         'ThEnter=%.3f ThExit=%.3f sig=%.3f jump=%d level=%d dip=%d stable=%d upL=%d upJ=%d dn=%d cool=%d hold=%d\n'], ...
        %         idnb, i, state_mode, fs_cur, Hacc_s(i), dHacc_raw(i), Th_enter, Th_exit, sigma_d, ...
        %         trig_up_jump, trig_up_level, big_dip, stable, up_count_level, up_count_jump, down_count, cooldown, hi_hold);
        % end
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
    logRec(idnb).dHacc = dHacc;
    logRec(idnb).BPM_est = BPM_est(1:frames);
    logRec(idnb).BPM0 = BPM0(1:frames);
    logRec(idnb).ACC_mag_25 = ACCmag25_log(1:frames);


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
selRecs = [6 15 17 21 22];
for idx = 1:numel(selRecs)
    r = selRecs(idx);
    if r <= numel(logRec) && ~isempty(logRec(r).BPM0)
        win_count = numel(logRec(r).Hacc);
        t_sec = (0:win_count-1) * step_sec; % window start times
        frames = 1:win_count;

        % HR comparison
        figure;
        plot(logRec(r).BPM0,'ro'); hold on; 
        plot(logRec(r).BPM_est,'o','Color','blue'); 
        xlabel('Time (frames)'); ylabel('HR (BPM)'); legend({'Ground truth','Estimates'});
        yyaxis right;
        plot(logRec(r).Hacc, '--', 'Color', 'green', 'LineWidth', 1.0, 'DisplayName', 'Hacc (25 Hz ACC)');
        plot(logRec(r).dHacc, '-', 'Color', [0 1 1], 'LineWidth', 1.0, 'DisplayName', 'dHacc (25 Hz ACC)');
        yline(0.15, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'DisplayName', 'Ref 0.15');
        yline(-0.15, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0, 'DisplayName', 'Ref -0.15');
        ylabel('Hacc / dHacc');
        yyaxis left;
        stairs(frames, logRec(r).FsUsed, '-'); ylabel('FsUsed (Hz)');
        title(sprintf('Recording %d (Adaptive Entropy Sampling)', r));
        grid on;

        % % Hacc and FsUsed over time (time axis in seconds, optional frame ticks on top)
        % figure;
        % yyaxis left; plot(frames, logRec(r).Hacc, '-'); ylabel('Hacc (bits)');
        % yyaxis right; stairs(frames, logRec(r).FsUsed, '-'); ylabel('FsUsed (Hz)');
        % xlabel('Time (frames)'); title(sprintf('Entropy & FsUsed - Recording %d', r));
        % grid on;

        % % Combined overlay: Hacc, dHacc_raw, FsUsed, ACC magnitude
        % figure;
        % ax1 = subplot(4,1,1); plot(frames, logRec(r).Hacc, '-'); ylabel('Hacc (bits)'); title(sprintf('Recording %d - Entropy & Sampling', r));
        % ax2 = subplot(4,1,2); plot(frames, logRec(r).dHacc, '-'); ylabel('dHacc (bits)'); grid on;
        % ax3 = subplot(4,1,3); stairs(frames, logRec(r).FsUsed, '-'); ylabel('FsUsed (Hz)'); grid on;
        % ax4 = subplot(4,1,4); plot(frames, logRec(r).ACC_mag_25, '-'); ylabel('ACC mag'); xlabel('Time (frames)'); grid on;
        % linkaxes([ax1 ax2 ax3 ax4],'x'); grid(ax1,'on');

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

function [BPM_val, state] = wfpv_one_frame_last(curDataRaw, fs0, fs_adc, fs_proc, FFTres, WFlength, searchHz, state, i, BPM_est, idnb, bpHz)
% Bandpass at 125 Hz, resample to fs_adc, then to 25 Hz internal, then run original WFPV logic
[b,a] = butter(4, bpHz/(fs0/2), 'bandpass');
curDataFilt = zeros(size(curDataRaw));
for c = 1:size(curDataRaw,1)
    curDataFilt(c,:) = filter(b,a,curDataRaw(c,:));
end
% resample to fs_adc
curData_adc = do_resample_last(curDataFilt, fs0, fs_adc);
% resample to internal 25 Hz if needed
if abs(fs_adc - fs_proc) < eps
    curData = curData_adc; fs = fs_proc;
else
    curData = do_resample_last(curData_adc, fs_adc, fs_proc);
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

function sig_res = do_resample_last(sig, fs_in, fs_out)
%DO_RESAMPLE Resample multi-channel signal (channels x samples) with proper anti-aliasing.
% Uses DECIMATE for integer downsampling ratios, otherwise RESAMPLE with rational p/q.
%
% sig:   [nCh x nSamp] (time along columns)
% fs_in: input sample rate (Hz)
% fs_out: output sample rate (Hz)

    if abs(fs_out - fs_in) < 1e-12
        sig_res = sig;
        return;
    end

    [nCh, nSamp] = size(sig);

    % ---------- Case 1: integer-factor downsample ----------
    r = fs_in / fs_out;
    if fs_out < fs_in && abs(r - round(r)) < 1e-9
        decim = round(r);
        % decimate includes an anti-alias lowpass; FIR is safer for bio signals
        outLens = ceil(nSamp / decim);
        sig_res = zeros(nCh, outLens);
        for c = 1:nCh
            y = decimate(sig(c,:), decim, 'fir');
            sig_res(c, 1:numel(y)) = y;
        end
        % trim in case decimate returned slightly longer
        sig_res = sig_res(:, 1:max(1, min(size(sig_res,2), outLens)));
        return;
    end

    % ---------- Case 2: integer-factor upsample ----------
    u = fs_out / fs_in;
    if fs_out > fs_in && abs(u - round(u)) < 1e-9
        up = round(u);
        % resample has a good interpolation/anti-imaging filter
        sig_res = resample(sig.', up, 1).';
        return;
    end

    % ---------- Case 3: rational resampling ----------
    ratio = fs_out / fs_in;
    [p,q] = rat(ratio, 1e-12);
    sig_res = resample(sig.', p, q).';
end
