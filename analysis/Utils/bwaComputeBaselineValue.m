function b = bwaComputeBaselineValue(vals, mask, transform)
%BWACOMPUTEBASELINEVALUE Compute baseline value for normalization.
% For log2_fold, this uses a geometric mean baseline (positive values only).

x = double(vals(mask));
x = x(isfinite(x));
if isempty(x)
    b = NaN;
    return
end

transform = lower(char(string(transform)));
if strcmp(transform, 'log2_fold')
    x = x(x > 0);
    if isempty(x)
        b = NaN;
        return
    end
    b = 2.^mean(log2(x), 'omitnan');
else
    b = mean(x, 'omitnan');
end
end
