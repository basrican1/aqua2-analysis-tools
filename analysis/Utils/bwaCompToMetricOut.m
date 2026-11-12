function mOut = bwaCompToMetricOut(comp, filesTbl, epochNames, titleName, normCfg)
%BWACOMPTOMETRICOUT Convert derived comp struct to canonical-like metric struct.

vals = double(comp.dffAll(:));
ep = double(comp.epMembIdx(:));
fid = strings(numel(vals),1);
cellID = strings(numel(vals),1);
if isfield(comp, 'fileID')
    ftmp = string(comp.fileID(:));
    fid(1:min(numel(fid), numel(ftmp))) = ftmp(1:min(numel(fid), numel(ftmp)));
end
if isfield(comp, 'memberlist')
    ctmp = string(comp.memberlist(:));
    cellID(1:min(numel(cellID), numel(ctmp))) = ctmp(1:min(numel(cellID), numel(ctmp)));
end
keep = isfinite(vals) & isfinite(ep) & ep > 0;
if any(fid ~= "")
    keep = keep & fid ~= "";
end
vals = vals(keep);
ep = ep(keep);
fid = fid(keep);
cellID = cellID(keep);
if isempty(vals)
    error('bwaCompToMetricOut:Empty', 'No valid rows in derived comp.');
end

if nargin >= 5 && isstruct(normCfg) && isfield(normCfg, 'mode') && ~strcmpi(string(normCfg.mode), "none")
    bIdx = 1;
    if isfield(normCfg, 'baselineEpochIdx') && ~isempty(normCfg.baselineEpochIdx)
        bIdx = double(normCfg.baselineEpochIdx);
    end
    vals = normalizeDerivedVals(vals, ep, cellID, normCfg.mode, normCfg.transform, bIdx);
end

T = table(fid, ep, vals, 'VariableNames', {'file_id','epoch_idx','value'});
T = groupsummary(T, {'file_id','epoch_idx'}, 'mean', 'value');
T.mean_value = T.mean_value;
T.median_value = T.mean_value;
T.std_value = zeros(height(T),1);
T.n_events = ones(height(T),1);
T.sem = zeros(height(T),1);
if isempty(epochNames)
    epochNamesCell = {};
else
    epochNamesCell = cellstr(string(epochNames(:)));
end
eNames = strings(height(T),1);
for i = 1:height(T)
    ei = T.epoch_idx(i);
    if ~isempty(epochNamesCell) && ei >= 1 && ei <= numel(epochNamesCell)
        eNames(i) = string(epochNamesCell{ei});
    else
        eNames(i) = string(['E' num2str(ei)]);
    end
end
T.epoch_name = eNames;

if ~isempty(filesTbl) && istable(filesTbl) && all(ismember({'file_id','animal_id'}, filesTbl.Properties.VariableNames))
    filesJ = filesTbl(:, {'file_id','animal_id'});
    filesJ.file_id = string(filesJ.file_id);
    T.file_id = string(T.file_id);
    J = outerjoin(T, filesJ, 'Keys', 'file_id', 'Type', 'left', 'MergeKeys', true);
    if ~ismember('animal_id', J.Properties.VariableNames)
        J.animal_id = repmat("", height(J), 1);
    end
else
    J = T;
    J.animal_id = J.file_id;
end

Te = table(fid, ep, vals, cellID, 'VariableNames', {'file_id','epoch_idx','value','primary_cell_uid'});
if ~isempty(filesTbl) && istable(filesTbl) && any(ismember({'file_id'}, filesTbl.Properties.VariableNames))
    cols = {'file_id'};
    for c = {'animal_id','slice_num'}
        if ismember(c{1}, filesTbl.Properties.VariableNames)
            cols{end+1} = c{1}; %#ok<AGROW>
        end
    end
    f2 = filesTbl(:, cols);
    f2.file_id = string(f2.file_id);
    Te.file_id = string(Te.file_id);
    Te = outerjoin(Te, f2, 'Keys', 'file_id', 'Type', 'left', 'MergeKeys', true);
end
if ~ismember('epoch_name', Te.Properties.VariableNames)
    eNamesEv = strings(height(Te),1);
    for i = 1:height(Te)
        ei = Te.epoch_idx(i);
        if ~isempty(epochNamesCell) && ei >= 1 && ei <= numel(epochNamesCell)
            eNamesEv(i) = string(epochNamesCell{ei});
        else
            eNamesEv(i) = string(['E' num2str(ei)]);
        end
    end
    Te.epoch_name = eNamesEv;
end

mOut = struct();
mOut.metricCol = 'value';
mOut.summaryByGroupEpoch = J;
mOut.eventValues = Te;
uEpoch = unique(J.epoch_idx, 'stable');
uNames = strings(numel(uEpoch),1);
for i = 1:numel(uEpoch)
    ei = uEpoch(i);
    if ~isempty(epochNamesCell) && ei >= 1 && ei <= numel(epochNamesCell)
        uNames(i) = string(epochNamesCell{ei});
    else
        uNames(i) = string(['E' num2str(ei)]);
    end
end
mOut.epochTable = table(uEpoch, uNames, 'VariableNames', {'epoch_idx','epoch_name'});
mOut.groupBy = {'animal_id','slice_num','file_id'};
mOut.titleName = titleName;
end

function out = normalizeDerivedVals(vals, ep, cellID, mode, transform, baseIdx)
out = NaN(size(vals));
mode = lower(char(string(mode)));
base = NaN(size(vals));

switch mode
    case 'by_cell_baseline_mean'
        u = unique(cellID, 'stable');
        for i = 1:numel(u)
            k = u(i);
            if strlength(k) == 0
                continue
            end
            idx = cellID == k;
            b = bwaComputeBaselineValue(vals, idx & ep == baseIdx, transform);
            base(idx) = b;
        end
    otherwise
        b = bwaComputeBaselineValue(vals, ep == baseIdx, transform);
        base(:) = b;
end

out = bwaApplyTransform(vals, base, transform);
end
