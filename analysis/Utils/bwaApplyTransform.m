function out = bwaApplyTransform(vals, base, transform)
%BWAAPPLYTRANSFORM Apply normalization transform to values given baseline.

out = NaN(size(vals));
vals = double(vals);
base = double(base);
transform = lower(char(string(transform)));

good = isfinite(vals) & isfinite(base);
switch transform
    case 'delta'
        out(good) = vals(good) - base(good);
    case 'percent_change'
        good = good & base ~= 0;
        out(good) = 100 * (vals(good) - base(good)) ./ base(good);
    otherwise % log2_fold
        good = good & vals > 0 & base > 0;
        out(good) = log2(vals(good) ./ base(good));
end
end
