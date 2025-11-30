%% ENTROPY_PROXY_LZ LZ-based entropy-density proxy (bits/symbol).
function H = entropy_proxy_lz(x, nbits)
    x = x(:)';
    if numel(x) < 8
        H = 0; return;
    end

    % Normalize + clip
    x = (x - mean(x)) / (std(x) + eps);
    x = max(min(x, 3), -3);

    % Quantize to alphabet A = 2^nbits over fixed [-3,3]
    A = 2^nbits;
    q = round((x + 3) * (A-1) / 6);
    q = min(max(q, 0), A-1);

    % Option A (recommended): use deltas so amplitude drifts don't dominate
    d = diff(q);
    % Map deltas to nonnegative symbols
    S = 2*(A-1) + 1;
    sym = d + (A-1);   % [0..S-1]
    n = numel(sym);

    if n < 8
        H = 0; return;
    end

    c = lz_phrase_count(sym);   % number of new phrases
    H = (c * log2(n)) / n;      % bits per symbol proxy
end


%% Helper: greedy "new phrase" count (shortest unseen substring parsing).
function c = lz_phrase_count(sym)
    sym = sym(:)';      % row
    n = numel(sym);
    c = 0;
    i = 1;

    % Use containers.Map for "seen phrases" (fine for n~200)
    dict = containers.Map('KeyType','char','ValueType','logical');

    while i <= n
        % Find smallest L such that sym(i:i+L-1) not seen before
        L = 1;
        while true
            j = i + L - 1;
            if j > n
                % leftover tail counts as a phrase
                c = c + 1;
                return;
            end
            key = phrase_key(sym(i:j));
            if ~isKey(dict, key)
                dict(key) = true;
                c = c + 1;
                i = j + 1;
                break;
            end
            L = L + 1;
        end
    end
end

function key = phrase_key(v)
    % Fast-ish stable serialization for small integers
    key = sprintf('%d,', v);
end