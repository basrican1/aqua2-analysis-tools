function labels = bwaMakeSliceLabels(fileIDs, filesTbl)
% Build readable, unique slice labels (preferred: yymmdd_S#) for file IDs.
%
% Inputs:
%   fileIDs  : vector/cellstr/string of file IDs.
%   filesTbl : optional canonical files table with file_id/exp_date/slice_num.

if nargin < 2
    filesTbl = table();
end
if istable(fileIDs)
    filesTbl = fileIDs;
    if ~ismember('file_id', filesTbl.Properties.VariableNames)
        error('bwaMakeSliceLabels:MissingFileID', ...
            'When the first input is a table, it must include a file_id column.');
    end
    fileIDs = filesTbl.file_id;
end

ids = string(fileIDs(:));
n = numel(ids);
labels = strings(n, 1);

for i = 1:n
    lbl = "";
    if istable(filesTbl) && ~isempty(filesTbl) && ...
            all(ismember({'file_id','exp_date','slice_num'}, filesTbl.Properties.VariableNames))
        idx = find(string(filesTbl.file_id) == ids(i), 1, 'first');
        if ~isempty(idx)
            d = regexprep(char(string(filesTbl.exp_date(idx))), '[^0-9]', '');
            if numel(d) >= 8
                d = d(end-5:end);
            end
            s = regexprep(char(string(filesTbl.slice_num(idx))), '[^0-9]', '');
            if isempty(s)
                s = char(string(filesTbl.slice_num(idx)));
            end
            if ~isempty(d) && ~isempty(s)
                lbl = string(sprintf('%s_S%s', d, s));
            end
        end
    end
    if strlength(lbl) == 0
        alt = string(bwaFormatFileIDLabel(char(ids(i))));
        if alt == ids(i) || strlength(alt) == 0
            alt = "f" + compose('%02d', i);
        end
        lbl = alt;
    end
    labels(i) = lbl;
end

labels = makeUnique(labels);
labels = cellstr(labels)';
end

function out = makeUnique(in)
out = in(:);
for i = 1:numel(out)
    if any(out(1:i-1) == out(i))
        k = 2;
        cand = out(i) + "_" + string(k);
        while any(out(1:i-1) == cand)
            k = k + 1;
            cand = out(i) + "_" + string(k);
        end
        out(i) = cand;
    end
end
end
