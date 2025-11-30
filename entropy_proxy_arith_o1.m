%% ENTROPY_PROXY_ARITH_O1  Order-1 adaptive model + ideal arithmetic length (bits/symbol)
function H = entropy_proxy_arith_o1(x, nbits)
    x = x(:)'; 
    if numel(x) < 4, H = 0; return; end

    % ---- Normalize + fixed-range quantization (stable across windows)
    x = (x - mean(x)) / (std(x) + eps);
    x = max(min(x, 3), -3);

    L = 2^nbits;
    q = round((x + 3) * (L-1) / 6);
    q = min(max(q, 0), L-1);

    % ---- Delta symbols -> alphabet [0..S-1]
    d = diff(q);
    S = 2*(L-1) + 1;
    sym = d + (L-1);        % [0..S-1]
    n = numel(sym);
    if n < 3, H = 0; return; end

    % ---- SCL-ish adaptive order-1 freq model
    % initialize with all ones (uniform prior) (like np.ones in your SCL snippet)
    counts = ones(S, S, 'double');  % counts(prev, curr)
    max_allowed_total = 2048;       % tune; prevents runaway totals like SCL
                                   % (you can set 512/1024/2048)
    bits = 0;
    prev = sym(1) + 1;             % MATLAB 1-index

    for t = 2:n
        curr = sym(t) + 1;

        row = counts(prev, :);
        p = row(curr) / sum(row);  % already has prior mass from ones()
        bits = bits - log2(p);

        % update after coding (important for decoder sync in real arithmetic coding)
        counts(prev, curr) = counts(prev, curr) + 1;

        % rescale row if needed (mirrors SCL’s “divide by 2, keep min 1” idea)
        if sum(counts(prev, :)) >= max_allowed_total
            counts(prev, :) = max(floor(counts(prev, :) / 2), 1);
        end

        prev = curr;
    end

    H = bits / (n-1);  % bits per delta-symbol
end