%% ENTROPY_PROXY_HIST 0th-order (histogram) entropy on delta symbols (bits/symbol).
function H = entropy_proxy_hist(x, nbits)
    x = x(:)';
    if numel(x) < 4
        H = 0; return;
    end

    % Normalize + clip to fixed range so entropy is comparable across windows
    x = (x - mean(x)) / (std(x) + eps);
    x = max(min(x, 3), -3);

    % Quantize to L levels over fixed [-3,3]
    L = 2^nbits;
    q = round((x + 3) * (L-1) / 6);
    q = min(max(q, 0), L-1);

    % Delta symbols (range roughly [-(L-1), +(L-1)])
    d = diff(q);
    S = 2*(L-1) + 1;            % alphabet size for deltas
    sym = d + (L-1);            % map to [0..S-1]

    % Histogram / plug-in entropy
    counts = accumarray(sym(:)+1, 1, [S 1]);
    p = counts / sum(counts);
    p = p(p>0);
    H = -sum(p .* log2(p));     % bits per delta-symbol
end