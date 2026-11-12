function out = mergeBWAPlotStyle(baseStyle, overrideStyle)
% Merge two style structs recursively.

if nargin < 1 || isempty(baseStyle)
    baseStyle = struct();
end
if nargin < 2 || isempty(overrideStyle)
    out = baseStyle;
    return
end
if ~isstruct(overrideStyle)
    error('mergeBWAPlotStyle:BadOverride', 'overrideStyle must be a struct.');
end

out = baseStyle;
f = fieldnames(overrideStyle);
for i = 1:numel(f)
    k = f{i};
    v = overrideStyle.(k);
    if isfield(out, k) && isstruct(out.(k)) && isstruct(v)
        out.(k) = mergeBWAPlotStyle(out.(k), v);
    else
        out.(k) = v;
    end
end
end
