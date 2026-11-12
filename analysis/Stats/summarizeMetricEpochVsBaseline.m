function out = summarizeMetricEpochVsBaseline(metricOut, varargin)
% Compute epoch-vs-baseline significance from a canonical-like metric output.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'metricOut', @isstruct);
addParameter(p, 'comparisonLevel', 'by_slice', @(x) ischar(x) || isstring(x));
addParameter(p, 'cellInclusion', 'complete_cells_only', @(x) ischar(x) || isstring(x));
addParameter(p, 'groupLabel', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'traceLabel', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'stat', 'mean', @(x) ischar(x) || isstring(x));
addParameter(p, 'sigAlpha', 0.05, @isnumeric);
addParameter(p, 'sigMethod', 'auto', @(x) ischar(x) || isstring(x));
parse(p, metricOut, varargin{:});

metricCol = metricOut.metricCol;
stat = lower(char(string(p.Results.stat)));
if strcmp(stat, 'median')
    ycol = ['median_' metricCol];
else
    ycol = ['mean_' metricCol];
end

[S, gcol, tvals, tname, ycol] = buildPlotTable(metricOut, lower(char(string(p.Results.comparisonLevel))), ...
    p.Results.groupLabel, p.Results.traceLabel, p.Results.cellInclusion, ycol);
TM = buildTraceEpochTable(S, tvals, gcol, ycol);
baselineIdx = min(double(metricOut.epochTable.epoch_idx));
sig = bwaComputeEpochVsBaselineSig(TM, baselineIdx, char(string(p.Results.sigMethod)), p.Results.sigAlpha);

if isfield(metricOut, 'epochTable') && istable(metricOut.epochTable) && ~isempty(metricOut.epochTable)
    sig = outerjoin(sig, metricOut.epochTable(:, {'epoch_idx','epoch_name'}), ...
        'Keys', 'epoch_idx', 'Type', 'left', 'MergeKeys', true);
end

out = struct();
out.metricCol = metricCol;
out.valueVar = ycol;
out.groupVar = gcol;
out.traceMode = tname;
out.summaryRows = S;
out.traceTable = TM;
out.significance = sig;
end

function [S, gcol, tvals, tname, ycol] = buildPlotTable(metricOut, comparisonLevel, groupLabel, traceLabel, cellInclusion, ycol)
S = metricOut.summaryByGroupEpoch;
if isempty(S)
    gcol = 'file_id';
    tvals = strings(0,1);
    tname = 'none';
    return
end

if isempty(groupLabel)
    gcol = 'animal_id';
else
    gcol = char(string(groupLabel));
end
if ~ismember(gcol, S.Properties.VariableNames)
    if ismember('animal_id', S.Properties.VariableNames)
        gcol = 'animal_id';
    else
        gcol = 'file_id';
    end
end

switch comparisonLevel
    case 'overall'
        tname = 'none';
        tvals = repmat("all", height(S), 1);
    case 'by_animal'
        tname = 'animal_id';
        if ismember('animal_id', S.Properties.VariableNames)
            tvals = string(S.animal_id);
            gcol = 'animal_id';
        else
            tvals = string(S.(gcol));
            tname = gcol;
        end
    case 'by_slice'
        if ismember('animal_id', S.Properties.VariableNames) && ismember('slice_num', S.Properties.VariableNames)
            tname = 'animal_slice';
            tvals = string(S.animal_id) + "|" + string(S.slice_num);
            if isempty(groupLabel)
                gcol = 'animal_id';
            end
        elseif ismember('file_id', S.Properties.VariableNames)
            tname = 'file_id';
            tvals = string(S.file_id);
        else
            tname = 'none';
            tvals = repmat("all", height(S), 1);
        end
    case 'by_cell'
        [S, gcol, ycol] = metricByCellSummary(metricOut, gcol, cellInclusion);
        if isempty(S)
            tname = 'none';
            tvals = strings(0,1);
            return
        end
        tname = 'primary_cell_uid';
        tvals = string(S.primary_cell_uid);
    otherwise
        req = char(string(traceLabel));
        if strcmpi(req, 'none') || isempty(req)
            tname = 'none';
            tvals = repmat("all", height(S), 1);
        elseif strcmpi(req, 'auto')
            if ismember('file_id', S.Properties.VariableNames)
                tname = 'file_id';
                tvals = string(S.file_id);
            else
                tname = gcol;
                tvals = string(S.(gcol));
            end
        elseif ismember(req, S.Properties.VariableNames)
            tname = req;
            tvals = string(S.(req));
        else
            tname = 'none';
            tvals = repmat("all", height(S), 1);
        end
end
end

function [S, gcol, ycol] = metricByCellSummary(metricOut, gcol, cellInclusion)
S = table();
ycol = ['mean_' metricOut.metricCol];
if ~isfield(metricOut, 'eventValues') || isempty(metricOut.eventValues)
    return
end
T = metricOut.eventValues;
if ~ismember('primary_cell_uid', T.Properties.VariableNames)
    if ismember('cell_uid', T.Properties.VariableNames)
        T.primary_cell_uid = string(T.cell_uid);
    else
        return
    end
end
mcol = metricOut.metricCol;
keep = isfinite(double(T.(mcol))) & T.epoch_idx > 0 & string(T.primary_cell_uid) ~= "";
T = T(keep, :);
if isempty(T)
    return
end
if ~ismember(gcol, T.Properties.VariableNames)
    if ismember('animal_id', T.Properties.VariableNames)
        gcol = 'animal_id';
    else
        gcol = 'file_id';
    end
end
S = groupsummary(T, {gcol,'primary_cell_uid','epoch_idx','epoch_name'}, 'mean', mcol);
S.(ycol) = S.(['mean_' mcol]);
S = expandCellEpochGrid(S, gcol, ycol, metricOut.epochTable);
S = applyCellInclusionPolicy(S, gcol, 'primary_cell_uid', ycol, metricOut.epochTable, cellInclusion);
S.sem = zeros(height(S),1);
end

function S = expandCellEpochGrid(Sin, gcol, ycol, epochTable)
if isempty(Sin) || isempty(epochTable)
    S = Sin;
    return
end
eIdx = epochTable.epoch_idx(:);
eName = string(epochTable.epoch_name(:));
grp = unique(Sin(:, {gcol, 'primary_cell_uid'}), 'rows', 'stable');
S = table();
for i = 1:height(grp)
    G = repmat(grp(i,:), numel(eIdx), 1);
    G.epoch_idx = eIdx;
    G.epoch_name = eName;
    J = outerjoin(G, Sin(:, {gcol,'primary_cell_uid','epoch_idx',ycol}), ...
        'Keys', {gcol,'primary_cell_uid','epoch_idx'}, 'Type', 'left', 'MergeKeys', true);
    if ~ismember('epoch_name', J.Properties.VariableNames)
        J.epoch_name = eName;
    end
    S = [S; J(:, {gcol,'primary_cell_uid','epoch_idx','epoch_name',ycol})]; %#ok<AGROW>
end
end

function S = applyCellInclusionPolicy(S, gcol, cellCol, ycol, epochTable, cellInclusion)
mode = lower(char(string(cellInclusion)));
if isempty(S)
    return
end
switch mode
    case 'available_case'
        return
    otherwise
        nEpoch = numel(unique(epochTable.epoch_idx));
        keys = unique(S(:, {gcol, cellCol}), 'rows', 'stable');
        keep = false(height(S),1);
        for i = 1:height(keys)
            idx = string(S.(gcol)) == string(keys.(gcol)(i)) & string(S.(cellCol)) == string(keys.(cellCol)(i));
            nPresent = sum(isfinite(double(S.(ycol)(idx))));
            if nPresent == nEpoch
                keep(idx) = true;
            end
        end
        S = S(keep, :);
end
end

function [x, y] = aggregateByEpoch(epochIdx, vals)
epochIdx = double(epochIdx(:));
vals = double(vals(:));
keep = isfinite(epochIdx) & isfinite(vals);
epochIdx = epochIdx(keep);
vals = vals(keep);
if isempty(epochIdx)
    x = zeros(0,1);
    y = zeros(0,1);
    return
end
[x, ~, gidx] = unique(epochIdx, 'stable');
y = accumarray(double(gidx), vals, [], @mean, NaN);
[x, ord] = sort(x);
y = y(ord);
end

function TM = buildTraceEpochTable(S, tvals, gcol, ycol)
tr = unique(tvals, 'stable');
TM = table(strings(0,1), strings(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', {'trace','group','epoch_idx','value'});
for i = 1:numel(tr)
    R = S(tvals == tr(i), :);
    if isempty(R)
        continue
    end
    [x, y] = aggregateByEpoch(R.epoch_idx, R.(ycol));
    t = repmat(tr(i), numel(x), 1);
    g = repmat(string(R.(gcol)(1)), numel(x), 1);
    TM = [TM; table(t, g, x, y, 'VariableNames', TM.Properties.VariableNames)]; %#ok<AGROW>
end
end
