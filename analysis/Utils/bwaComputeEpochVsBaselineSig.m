function sig = bwaComputeEpochVsBaselineSig(TM, baselineIdx, method, alpha)
% Compute epoch-vs-baseline significance from trace-level epoch summaries.
%
% TM must contain:
%   trace, epoch_idx, value

uEpoch = unique(double(TM.epoch_idx), 'stable');
sig = table(uEpoch, nan(numel(uEpoch),1), false(numel(uEpoch),1), ...
    'VariableNames', {'epoch_idx','p_value','is_sig'});

for i = 1:numel(uEpoch)
    e = uEpoch(i);
    if e == baselineIdx
        continue
    end
    B = TM(TM.epoch_idx == baselineIdx, :);
    E = TM(TM.epoch_idx == e, :);
    if isempty(B) || isempty(E)
        continue
    end
    [x, y] = matchTraceValues(B, E);
    if numel(x) >= 3 && numel(y) >= 3
        p = doSigTest(x, y, method);
    else
        p = doSigTest(double(B.value), double(E.value), 'unpaired');
    end
    sig.p_value(i) = p;
    sig.is_sig(i) = isfinite(p) && p < alpha;
end
end

function [x, y] = matchTraceValues(B, E)
tb = table(string(B.trace), double(B.value), 'VariableNames', {'trace','val_b'});
te = table(string(E.trace), double(E.value), 'VariableNames', {'trace','val_e'});
M = innerjoin(tb, te, 'Keys', 'trace');
x = double(M.val_b);
y = double(M.val_e);
keep = isfinite(x) & isfinite(y);
x = x(keep);
y = y(keep);
end

function p = doSigTest(x, y, method)
x = double(x(:));
y = double(y(:));
x = x(isfinite(x));
y = y(isfinite(y));
p = NaN;
if isempty(x) || isempty(y)
    return
end

m = lower(char(string(method)));
if strcmp(m, 'auto')
    m = 'paired';
end

switch m
    case 'paired'
        n = min(numel(x), numel(y));
        x = x(1:n);
        y = y(1:n);
        if n < 3
            return
        end
        d = y - x;
        try
            p = signrank(d, 0, 'method', 'approximate');
        catch
            try
                [~, p] = ttest(d, 0);
            catch
                p = NaN;
            end
        end
    otherwise
        if numel(x) < 2 || numel(y) < 2
            return
        end
        try
            [~, p] = ttest2(x, y);
        catch
            p = ranksum(x, y, 'method', 'approximate');
        end
end
end
