function opts = bwa_aqua2_read_preset(csvPath)
% Read an AQuA2 GUI option export without using its recorded input path.

raw = readcell(csvPath, 'Delimiter', ',');
if size(raw, 1) < 2
    error('AQuA2:InvalidPreset', 'Preset must contain a header and value row: %s', csvPath);
end

names = string(raw(1, :));
values = raw(2, :);
opts = struct();

% Columns 1:10 describe the old input and its size. The current input is
% loaded separately, so only named algorithm settings are imported.
for idx = 11:numel(names)
    name = strtrim(names(idx));
    if name == "" || ismissing(name)
        continue
    end
    value = values{idx};
    if ismissingValue(value)
        continue
    end
    if ischar(value) || isstring(value)
        numericValue = str2double(value);
        if ~isnan(numericValue) || strcmpi(strtrim(string(value)), "NaN")
            value = numericValue;
        end
    end
    opts.(char(name)) = value;
end
end

function tf = ismissingValue(value)
tf = isempty(value);
if ismissing(value)
    tf = true;
elseif ischar(value) || isstring(value)
    tf = strlength(strtrim(string(value))) == 0;
end
end
