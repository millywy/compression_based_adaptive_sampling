function [FsNext, st] = ctrl_entropy_online(Hacc_bits, accMag, st)
    % --- init
    if isempty(st)
        st.mode = 6.25;                 % current mode (Hz)
        st.dwell = 0;                 % windows spent in current mode
        st.cool  = 0;                 % cooldown windows remaining
        st.muH = 0.5; st.madH = 0.05; % online baseline for Hnorm
        st.lowCount = 0;              % consecutive "safe to go down"
    end

    % --- normalize entropy (nbits=2 => S=7 deltas)
    S = 7;
    Hn = Hacc_bits / log2(S);         % ~[0,1]

    % --- update online baseline (slow)
    beta = 0.02;                      % ~50-window time constant
    st.muH  = (1-beta)*st.muH  + beta*Hn;
    st.madH = (1-beta)*st.madH + beta*abs(Hn - st.muH);
    zH = (Hn - st.muH) / (st.madH + 1e-6);

    % --- smooth the decision variable a bit (optional)
    if ~isfield(st,'zHsm'), st.zHsm = zH; end
    st.zHsm = 0.7*st.zHsm + 0.3*zH;

    % --- dwell/cooldown bookkeeping
    st.dwell = st.dwell + 1;
    if st.cool > 0
        st.cool = st.cool - 1;
        FsNext = st.mode;
        return;
    end

    % --- thresholds in "relative units" (no percentiles needed)
    % upshift if entropy is much worse than baseline OR accMag is high
    up = (st.zHsm > 1.0) || (accMag > 0.8);

    % downshift candidate if entropy is better than baseline AND accMag low-ish
    down_ok = (st.zHsm < -0.7) && (accMag < 0.5);

    % --- min dwell to prevent thrash
    minDwell = 4;    % hold each state at least 4 windows (32s if 8s windows)

    if up && (st.mode < 25) && (st.dwell >= 1)   % allow quick upshift
        st.mode = 25;
        st.dwell = 0;
        st.cool = 2;                              % ignore switches for 2 windows
        st.lowCount = 0;

    else
        if down_ok
            st.lowCount = st.lowCount + 1;
        else
            st.lowCount = 0;
        end

        if (st.lowCount >= 4) && (st.dwell >= minDwell)
            % step down one level only
            if st.mode == 25
                st.mode = 12;   % MID
            elseif st.mode == 12
                st.mode = 6;    % LOW
            end
            st.dwell = 0;
            st.cool = 2;
            st.lowCount = 0;
        end
    end

    FsNext = st.mode;
end