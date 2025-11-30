function [FsNext, st] = ctrl_entropy_2state(Hacc_bits, accMag, nbits_entropy, st)
%CTRL_ENTROPY_2STATE Online 2-state FS controller (25 Hz <-> 12.5 Hz)
% Inputs:
%   Hacc_bits     - entropy proxy (bits) computed from ACC features per window
%   accMag        - ACC magnitude (per window), already in your plots
%   nbits_entropy - the bit depth used in your entropy quantizer (e.g., 2)
%   st            - persistent state struct (pass [] on first call)
% Output:
%   FsNext        - {25, 12.5} for the NEXT window
%   st            - updated state

    % ====== Params you may tune (but these defaults are sane) ======
    FS_HI = 25;
    FS_LO = 12.5;

    % Dwell/cooldown counts in windows (8s window -> 4 = 32s)
    MIN_DWELL_HI = 2;    % keep HI at least this long once you switch up
    MIN_DWELL_LO = 4;    % keep LO at least this long once you switch down
    COOLDOWN     = 2;    % after any switch, ignore switches for this many windows

    % Require sustained evidence to go down
    DOWN_M        = 3;   % consecutive "safe" windows before dropping to LO

    % Online baseline update speed for entropy (fully online; no percentiles needed)
    BETA = 0.02;         % ~50-window time constant

    % Relative thresholds in z-score units
    Z_UP   =  +1.0;      % if entropy worse than baseline by this much -> go HI
    Z_DOWN =  -0.8;      % if entropy better than baseline by this much -> consider LO

    % ACC magnitude safety thresholds (you can tune these to your accMag scaling)
    A_UP   = 0.85;       % strong motion -> HI
    A_DOWN = 0.55;       % only allow LO when motion is below this

    % ====== Init state ======
    if nargin < 4 || isempty(st)
        st.mode = FS_HI;           % start high (safe)
        st.dwell = 0;
        st.cool = 0;
        st.muH = 0.5;
        st.madH = 0.05;
        st.zHsm = 0;
        st.downCount = 0;
        st.reacq = 0;              % optional: reacquisition burst counter
    end

    % ====== Normalize entropy to [0,1] range-ish ======
    % delta alphabet size: S = 2*(2^nbits-1) + 1
    L = 2^nbits_entropy;
    S = 2*(L-1) + 1;
    Hn = Hacc_bits / log2(S);      % normalized entropy (dimensionless)

    % ====== Online baseline (mu + MAD) ======
    st.muH  = (1-BETA)*st.muH  + BETA*Hn;
    st.madH = (1-BETA)*st.madH + BETA*abs(Hn - st.muH);
    zH = (Hn - st.muH) / (st.madH + 1e-6);

    % Smooth z-score to reduce flicker
    st.zHsm = 0.7*st.zHsm + 0.3*zH;

    % ====== Bookkeeping ======
    st.dwell = st.dwell + 1;
    if st.cool > 0
        st.cool = st.cool - 1;
        FsNext = st.mode;
        return;
    end

    % Optional: reacquisition burst after going HI (keeps estimator stable)
    if st.reacq > 0
        st.reacq = st.reacq - 1;
        st.mode = FS_HI;
        FsNext = st.mode;
        return;
    end

    % ====== Decision logic ======
    up_condition = (st.zHsm > Z_UP) || (accMag > A_UP);
    down_ok      = (st.zHsm < Z_DOWN) && (accMag < A_DOWN);

    if st.mode == FS_LO
        % In LOW: go HIGH quickly if things look bad, but respect min dwell a bit
        if up_condition && st.dwell >= MIN_DWELL_LO
            st.mode  = FS_HI;
            st.dwell = 0;
            st.cool  = COOLDOWN;
            st.downCount = 0;
            st.reacq = 1;      % 1 extra HI window for reacq after switching up
        end

    else
        % In HIGH: only go LOW if down condition holds for DOWN_M windows AND min dwell met
        if down_ok
            st.downCount = st.downCount + 1;
        else
            st.downCount = 0;
        end

        if (st.downCount >= DOWN_M) && (st.dwell >= MIN_DWELL_HI)
            st.mode  = FS_LO;
            st.dwell = 0;
            st.cool  = COOLDOWN;
            st.downCount = 0;
        end
    end

    FsNext = st.mode;
end