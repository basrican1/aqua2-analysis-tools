function out = bwaFormatFileIDLabel(fileID)
% Format file IDs into compact slice labels: yymmdd_S# when possible.
% Falls back to original fileID text if tokens cannot be parsed.

s = char(string(fileID));
s = strtrim(s);
if isempty(s)
    out = s;
    return
end

% Date token: contiguous yyyymmdd / yymmdd, or separated yyyy-mm-dd / yyyy_mm_dd.
dateTok = regexp(s, '(20\d{6}|\d{6})', 'match', 'once');
if isempty(dateTok)
    d = regexp(s, '(20\d{2})[\-_\s]?(\d{2})[\-_\s]?(\d{2})', 'tokens', 'once');
    if ~isempty(d)
        dateTok = [d{1} d{2} d{3}];
    end
end

% Slice token: "slice1", "sl1", "s1", or segment token like "_S1_".
sliceTok = regexp(s, '(?:slice|sl|s)[_\-\s]*(\d+)', 'tokens', 'once', 'ignorecase');
if isempty(sliceTok)
    sliceTok = regexp(s, '(?:^|[_\-\s])s(\d+)(?:$|[_\-\s])', 'tokens', 'once', 'ignorecase');
end

if isempty(dateTok) || isempty(sliceTok)
    out = s;
    return
end

if numel(dateTok) >= 8
    datePart = dateTok(end-5:end);
elseif numel(dateTok) == 6
    datePart = dateTok;
else
    out = s;
    return
end

sNum = str2double(sliceTok{1});
if ~isfinite(sNum) || isnan(sNum)
    out = s;
    return
end

out = sprintf('%s_S%d', datePart, sNum);
end
