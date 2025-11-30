% Inputs per window: Hacc (bits), accMag (scalar)
% Output: FsUsed for NEXT window

function [FsNext, state] = entropy_controller(Hacc, accMag, state)

    % ---- normalize entropy
    % For nbits=2: L=4, S=7  (in general: S = 2*(2^nbits - 1) + 1)
    S = 7;
    Hnorm = Hacc / log2(S);

    % ---- thresholds (start here, then tune by percentiles)
    H_low  = 0.45;   H_high = 0.58;
    A_low  = 0.35;   A_high = 0.75;

    % ---- dwell params
    M_down = 3;
    if ~isfield(state,'downCount'), state.downCount = 0; end
    if ~isfield(state,'mode'), state.mode = "MID"; end  % "LOW","MID","HIGH"

    % ---- upshift (immediate, OR)
    if (Hnorm > H_high) || (accMag > A_high)
        state.mode = "HIGH";
        state.downCount = 0;

    else
        % ---- downshift candidate (AND, sustained)
        if (Hnorm < H_low) && (accMag < A_low)
            state.downCount = state.downCount + 1;
        else
            state.downCount = 0;
        end

        if state.downCount >= M_down
            state.mode = "LOW";
        else
            state.mode = "MID";
        end
    end

    % ---- map mode -> Fs
    switch state.mode
        case "LOW",  FsNext = 6.25;
        case "MID",  FsNext = 12.5;
        case "HIGH", FsNext = 25;
    end
end