function dataType = datTreeToDataType(datTree)
switch char(string(datTree))
    case 'basic.area'
        dataType = 'area';
    case 'curve.duration'
        dataType = 'duration';
    case 'curve.rise19'
        dataType = 'rise';
    case 'curve.fall91'
        dataType = 'fall';
    case 'curve.decayTau'
        dataType = 'decayTau';
    case 'curve.dffAUC'
        dataType = 'AUC';
    case 'curve.dffMax2'
        dataType = 'amplitude';
    case 'propagation.maxPropSpeed'
        dataType = 'maxspeed';
    otherwise
        dataType = '';
end
end

function BWA_sub = reduceBWAByFiles(BWA, fileIDs, filesTbl, eventsTbl)
if isempty(fileIDs)
    BWA_sub = BWA;
    return
end
ids = cellstr(string(fileIDs(:)));
BWA_sub = ReduceAQ2ByFile(BWA, ids, 'filesTbl', filesTbl, 'eventsTbl', eventsTbl);
end

function memberlist = computeMemberlistForBWA(BWA, opts, DB)
[memberlist, ~, ~, ~] = eventBelongsCell(BWA.fts1, opts, DB);
end

function polSub = subsetPolarsData(polData, fileIDs)
mask = ismember(polData.uniqueFiles(:), fileIDs(:));
if ~any(mask)
    error('subsetPolarsData:NoMatch', 'No matching file IDs found in polar data.');
end
polSub = struct();
fields = fieldnames(polData);
for i = 1:numel(fields)
    f = fields{i};
    val = polData.(f);
    if ischar(val) || isstring(val)
        polSub.(f) = val;
        continue
    end
    if iscell(val) && numel(val) == numel(polData.uniqueFiles)
        polSub.(f) = val(mask);
        continue
    end
    if isnumeric(val) && size(val, 1) == numel(polData.uniqueFiles)
        polSub.(f) = val(mask, :);
        continue
    end
    polSub.(f) = val;
end
polSub.uniqueFiles = polData.uniqueFiles(mask);
polSub.subplotGrid = [1 numel(polSub.uniqueFiles)];
end
