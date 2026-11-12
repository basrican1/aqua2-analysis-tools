function polData = augment_polars_with_cell_metrics(result, polData)
% Append selected table-native cell and slice metrics into spider payload.
%
% Uses result.metrics.registry to discover selected local metrics and appends
% per-file ExperimentalCondition-summary values to polData.valuesA / valuesC and labels.

if ~isstruct(polData) || ~all(isfield(polData, {'labels','valuesA','valuesC','uniqueFiles'}))
    return
end
if ~isstruct(result) || ~isfield(result, 'metrics') || ~isstruct(result.metrics) || ...
        ~isfield(result.metrics, 'registry') || ~istable(result.metrics.registry) || ...
        ~isfield(result.metrics, 'tables') || ~isstruct(result.metrics.tables)
    return
end

reg = result.metrics.registry;
files = string(polData.uniqueFiles(:));
if isempty(files)
    return
end

tables = result.metrics.tables;
labels = cellstr(string(polData.labels(:)));
A = polData.valuesA;
C = polData.valuesC;
entityA = getEntityValueCells(polData, 'entityValuesA', numel(labels));
entityC = getEntityValueCells(polData, 'entityValuesC', numel(labels));
entityAByFile = getEntityValueCells(polData, 'entityValuesAByFile', numel(labels));
entityCByFile = getEntityValueCells(polData, 'entityValuesCByFile', numel(labels));
allowedCell = defaultSpiderCellMetricIDs();
allowedSlice = defaultSpiderSliceMetricIDs();

for i = 1:height(reg)
    metricID = string(reg.metric_id(i));
    if ~(startsWith(metricID, "cell.") || startsWith(metricID, "slice."))
        continue
    end
    if metricID == "cell.frequency"
        % Existing polars already includes frequency branch.
        continue
    end
    if startsWith(metricID, "cell.") && ~any(metricID == allowedCell)
        continue
    end
    if startsWith(metricID, "slice.") && ~any(metricID == allowedSlice)
        continue
    end

    cVals = computeRegistryBranchValues(metricID, reg(i,:), tables, files);
    if isempty(cVals)
        continue
    end
    cEntity = computeRegistryBranchEntityValues(metricID, reg(i,:), tables);
    cEntityByFile = computeRegistryBranchEntityValuesByFile(metricID, reg(i,:), tables, files);

    if startsWith(metricID, "cell.")
        aVals = computeEventAnalogBranchValues(metricID, tables, files);
        aEntity = computeEventAnalogEntityValues(metricID, tables);
        aEntityByFile = computeEventAnalogEntityValuesByFile(metricID, tables, files);
        if isempty(aVals)
            aVals = cVals;
        end
        if isempty(aEntity)
            aEntity = cEntity;
        end
        if isempty(aEntityByFile)
            aEntityByFile = cEntityByFile;
        end
    else
        aVals = nan(size(cVals));
        aEntity = [];
        aEntityByFile = cell(numel(files), 1);
    end

    metricLabel = char(strrep(metricID, '.', '_'));
    if any(strcmp(metricLabel, labels))
        continue
    end

    labels{end+1,1} = metricLabel; %#ok<AGROW>
    A(:, end+1) = aVals; %#ok<AGROW>
    C(:, end+1) = cVals; %#ok<AGROW>
    entityA{end+1,1} = aEntity; %#ok<AGROW>
    entityC{end+1,1} = cEntity; %#ok<AGROW>
    entityAByFile{end+1,1} = aEntityByFile; %#ok<AGROW>
    entityCByFile{end+1,1} = cEntityByFile; %#ok<AGROW>
    if isfield(polData, 'sigMaskA') && isnumeric(polData.sigMaskA) && size(polData.sigMaskA, 1) == size(A, 1) && size(polData.sigMaskA, 2) == numel(labels)-1
        polData.sigMaskA(:, end+1) = false; %#ok<AGROW>
    end
    if isfield(polData, 'sigMaskC') && isnumeric(polData.sigMaskC) && size(polData.sigMaskC, 1) == size(C, 1) && size(polData.sigMaskC, 2) == numel(labels)-1
        polData.sigMaskC(:, end+1) = false; %#ok<AGROW>
    end
    if isfield(polData, 'panelLimitsLog2') && isnumeric(polData.panelLimitsLog2) && numel(polData.panelLimitsLog2) == numel(labels)-1
        polData.panelLimitsLog2(end+1) = NaN; %#ok<AGROW>
    end
end

if isfield(polData, 'sigMaskA') && isnumeric(polData.sigMaskA) && size(polData.sigMaskA,2) < numel(labels)
    polData.sigMaskA(:, end+1:numel(labels)) = false;
end
if isfield(polData, 'sigMaskC') && isnumeric(polData.sigMaskC) && size(polData.sigMaskC,2) < numel(labels)
    polData.sigMaskC(:, end+1:numel(labels)) = false;
end
if isfield(polData, 'panelLimitsLog2') && isnumeric(polData.panelLimitsLog2) && numel(polData.panelLimitsLog2) < numel(labels)
    polData.panelLimitsLog2(end+1:numel(labels)) = NaN;
end

polData.labels = labels;
polData.valuesA = A;
polData.valuesC = C;
polData.entityValuesA = entityA;
polData.entityValuesC = entityC;
polData.entityValuesAByFile = entityAByFile;
polData.entityValuesCByFile = entityCByFile;
end


function vals = computeRegistryBranchValues(metricID, regRow, tables, files)
tname = char(string(regRow.source_table));
vals = [];
if ~isfield(tables, tname)
    return
end
T = tables.(tname);
if ~istable(T) || isempty(T) || ~all(ismember({'file_id','epoch_idx'}, T.Properties.VariableNames))
    return
end
normCol = char(string(regRow.value_column_norm));
rawCol = char(string(regRow.value_column_raw));
valCol = '';
if strlength(string(normCol)) > 0 && ismember(normCol, T.Properties.VariableNames)
    valCol = normCol;
elseif ismember(rawCol, T.Properties.VariableNames)
    valCol = rawCol;
end
if isempty(valCol)
    return
end
vals = nan(numel(files),1);
for k = 1:numel(files)
    F = T(string(T.file_id) == files(k), :);
    vals(k) = meanCnoPair(F, valCol);
end
end

function vals = computeRegistryBranchEntityValues(metricID, regRow, tables)
tname = char(string(regRow.source_table));
vals = [];
if startsWith(string(metricID), "slice.")
    return
end
if ~isfield(tables, tname)
    return
end
T = tables.(tname);
if ~istable(T) || isempty(T) || ~all(ismember({'epoch_idx'}, T.Properties.VariableNames))
    return
end
normCol = char(string(regRow.value_column_norm));
rawCol = char(string(regRow.value_column_raw));
valCol = '';
if strlength(string(normCol)) > 0 && ismember(normCol, T.Properties.VariableNames)
    valCol = normCol;
elseif ismember(rawCol, T.Properties.VariableNames)
    valCol = rawCol;
end
if isempty(valCol)
    return
end
vals = combinedPostEntityValues(T, valCol, chooseEntityKey(T, 'cell'));
end

function valsByFile = computeRegistryBranchEntityValuesByFile(metricID, regRow, tables, files)
tname = char(string(regRow.source_table));
valsByFile = cell(numel(files), 1);
if startsWith(string(metricID), "slice.")
    return
end
if ~isfield(tables, tname)
    return
end
T = tables.(tname);
if ~istable(T) || isempty(T) || ~all(ismember({'file_id','epoch_idx'}, T.Properties.VariableNames))
    return
end
normCol = char(string(regRow.value_column_norm));
rawCol = char(string(regRow.value_column_raw));
valCol = '';
if strlength(string(normCol)) > 0 && ismember(normCol, T.Properties.VariableNames)
    valCol = normCol;
elseif ismember(rawCol, T.Properties.VariableNames)
    valCol = rawCol;
end
if isempty(valCol)
    return
end
entityKey = chooseEntityKey(T, 'cell');
for k = 1:numel(files)
    F = T(string(T.file_id) == files(k), :);
    valsByFile{k} = combinedPostEntityValues(F, valCol, entityKey);
end
end

function vals = computeEventAnalogBranchValues(metricID, tables, files)
vals = [];
[tname, valCol, agg] = eventAnalogSpec(metricID, tables);
if isempty(tname) || ~isfield(tables, tname)
    return
end
if strlength(string(valCol)) == 0
    return
end
T = tables.(tname);
if ~istable(T) || isempty(T) || ~all(ismember({'file_id','epoch_idx'}, T.Properties.VariableNames))
    return
end
vals = nan(numel(files),1);
for k = 1:numel(files)
    F = T(string(T.file_id) == files(k), :);
    vals(k) = agg(F, valCol);
end
end

function vals = computeEventAnalogEntityValues(metricID, tables)
vals = [];
[tname, valCol, agg] = eventAnalogSpec(metricID, tables); %#ok<ASGLU>
if isempty(tname) || ~isfield(tables, tname)
    return
end
if strlength(string(valCol)) == 0
    return
end
T = tables.(tname);
if ~istable(T) || isempty(T) || ~ismember('epoch_idx', T.Properties.VariableNames)
    return
end
switch char(metricID)
    case 'cell.frac_events_with_neighbor_recruitment'
        vals = combinedPostValuesByEventAny(T, valCol);
    case 'cell.mean_neighbor_peak_overlap_frac_of_cell'
        vals = combinedPostValuesByEventMean(T, valCol, 'is_recruited_partition', true);
    case 'cell.mean_neighbor_partitions_per_event'
        vals = combinedPostValuesByEventCount(T, valCol);
    otherwise
        vals = combinedPostEntityValues(T, valCol, chooseEntityKey(T, 'event'));
end
end

function valsByFile = computeEventAnalogEntityValuesByFile(metricID, tables, files)
valsByFile = cell(numel(files), 1);
[tname, valCol, agg] = eventAnalogSpec(metricID, tables); %#ok<ASGLU>
if isempty(tname) || ~isfield(tables, tname)
    return
end
if strlength(string(valCol)) == 0
    return
end
T = tables.(tname);
if ~istable(T) || isempty(T) || ~all(ismember({'file_id','epoch_idx'}, T.Properties.VariableNames))
    return
end
for k = 1:numel(files)
    F = T(string(T.file_id) == files(k), :);
    switch char(metricID)
        case 'cell.frac_events_with_neighbor_recruitment'
            valsByFile{k} = combinedPostValuesByEventAny(F, valCol);
        case 'cell.mean_neighbor_peak_overlap_frac_of_cell'
            valsByFile{k} = combinedPostValuesByEventMean(F, valCol, 'is_recruited_partition', true);
        case 'cell.mean_neighbor_partitions_per_event'
            valsByFile{k} = combinedPostValuesByEventCount(F, valCol);
        otherwise
            valsByFile{k} = combinedPostEntityValues(F, valCol, chooseEntityKey(F, 'event'));
    end
end
end

function [tname, valCol, agg] = eventAnalogSpec(metricID, tables)
tname = '';
valCol = '';
agg = @meanCnoPair;
switch char(metricID)
    case 'cell.frac_multicell_events'
        tname = 'event_spread';
        valCol = pickFirstAvailableColumn(tables, tname, {'is_multicell_norm','is_multicell'});
    case 'cell.mean_cells_per_multicell_event'
        tname = 'event_spread';
        valCol = pickFirstAvailableColumn(tables, tname, {'spread_count_norm','spread_count'});
    case 'cell.frac_soma_events'
        tname = 'event_soma';
        valCol = pickFirstAvailableColumn(tables, tname, {'soma_covered_ge50_norm','soma_covered_ge50'});
    case 'cell.frac_toward_events'
        tname = 'event_soma';
        valCol = pickFirstAvailableColumn(tables, tname, {'is_toward_event_norm','is_toward_event'});
    case 'cell.frac_away_events'
        tname = 'event_soma';
        valCol = pickFirstAvailableColumn(tables, tname, {'is_away_event_norm','is_away_event'});
    case 'cell.frac_traverse_events'
        tname = 'event_soma';
        valCol = pickFirstAvailableColumn(tables, tname, {'is_traverse_event_norm','is_traverse_event'});
    case 'cell.mean_motion_index'
        tname = 'event_soma';
        valCol = pickFirstAvailableColumn(tables, tname, {'motion_index_norm','motion_index'});
    case 'cell.mean_soma_centroid_dist_um'
        tname = 'event_soma';
        valCol = pickFirstAvailableColumn(tables, tname, {'soma_centroid_dist_um_norm','soma_centroid_dist_um'});
    case 'cell.mean_soma_start_dist_um'
        tname = 'event_soma';
        valCol = pickFirstAvailableColumn(tables, tname, {'soma_start_dist_um_norm','soma_start_dist_um'});
    case 'cell.mean_soma_end_dist_um'
        tname = 'event_soma';
        valCol = pickFirstAvailableColumn(tables, tname, {'soma_end_dist_um_norm','soma_end_dist_um'});
    case 'cell.frac_net_toward_soma_events'
        tname = 'event_soma';
        valCol = pickFirstAvailableColumn(tables, tname, {'is_net_toward_soma_norm','is_net_toward_soma'});
    case 'cell.frac_net_away_soma_events'
        tname = 'event_soma';
        valCol = pickFirstAvailableColumn(tables, tname, {'is_net_away_soma_norm','is_net_away_soma'});
    case 'cell.mean_soma_net_distance_change_um'
        tname = 'event_soma';
        valCol = pickFirstAvailableColumn(tables, tname, {'soma_net_distance_change_um_norm','soma_net_distance_change_um'});
    case 'cell.mean_perimeter_um'
        tname = 'event_shape';
        valCol = pickFirstAvailableColumn(tables, tname, {'perimeter_um_norm','perimeter_um'});
    case 'cell.mean_circularity'
        tname = 'event_shape';
        valCol = pickFirstAvailableColumn(tables, tname, {'circularity_norm','circularity'});
    case 'cell.frac_propagative_events'
        tname = 'event_propagation';
        valCol = pickFirstAvailableColumn(tables, tname, {'is_propagative_norm','is_propagative'});
    case 'cell.mean_propagation_distance_um'
        tname = 'event_propagation';
        valCol = pickFirstAvailableColumn(tables, tname, {'propagation_distance_um_norm','propagation_distance_um'});
    case 'cell.mean_frac_outside_origin_cell'
        tname = 'event_boundary';
        valCol = pickFirstAvailableColumn(tables, tname, {'frac_outside_cell_norm','frac_outside_cell'});
    case 'cell.frac_escape_events_gt_threshold'
        tname = 'event_boundary';
        valCol = pickFirstAvailableColumn(tables, tname, {'is_escape_gt_threshold_norm','is_escape_gt_threshold'});
    case 'cell.frac_events_with_neighbor_recruitment'
        tname = 'event_partition';
        valCol = pickFirstAvailableColumn(tables, tname, {'is_recruited_partition_norm','is_recruited_partition'});
        agg = @meanCnoPairByEvent;
    case 'cell.mean_neighbor_entry_delay_sec'
        tname = 'event_partition';
        valCol = pickFirstAvailableColumn(tables, tname, {'cell_entry_delay_from_primary_sec_norm','cell_entry_delay_from_primary_sec'});
        agg = @meanCnoPair;
    case 'cell.mean_neighbor_peak_overlap_frac_of_cell'
        tname = 'event_partition';
        valCol = pickFirstAvailableColumn(tables, tname, {'cell_peak_overlap_frac_of_cell_norm','cell_peak_overlap_frac_of_cell'});
        agg = @meanCnoPairRecruitedOnly;
    case 'cell.mean_neighbor_partitions_per_event'
        tname = 'event_partition';
        valCol = pickFirstAvailableColumn(tables, tname, {'is_recruited_partition_norm','is_recruited_partition'});
        agg = @meanCnoPairRecruitedCountByEvent;
end
end

function name = pickFirstAvailableColumn(tables, tname, candidates)
name = '';
if ~isfield(tables, tname)
    return
end
T = tables.(tname);
for i = 1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        name = candidates{i};
        return
    end
end
end

function v = meanCnoPairByEvent(T, col)
if isempty(T) || ~ismember(col, T.Properties.VariableNames) || ~ismember('event_uid', T.Properties.VariableNames)
    v = NaN;
    return
end
vals = nan(1,3);
for e = 1:3
    Te = T(T.epoch_idx == e, {'event_uid', col});
    if isempty(Te)
        continue
    end
    [G,~] = findgroups(string(Te.event_uid));
    ev = splitapply(@(x) any(logical(x(isfinite(x)))), double(Te.(col)), G);
    ev = double(ev(:));
    ev = ev(isfinite(ev));
    if ~isempty(ev)
        vals(e) = mean(ev, 'omitnan');
    end
end
if all(isnan(vals(2:3)))
    v = NaN;
elseif isnan(vals(3))
    v = vals(2);
elseif isnan(vals(2))
    v = vals(3);
else
    v = mean(vals(2:3), 'omitnan');
end
end
function v = meanCnoPair(T, col)
if isempty(T) || ~ismember(col, T.Properties.VariableNames)
    v = NaN;
    return
end
vals = nan(1,3);
for e = 1:3
    x = double(T.(col)(T.epoch_idx == e));
    x = x(isfinite(x));
    if ~isempty(x)
        vals(e) = mean(x, 'omitnan');
    end
end
if all(isnan(vals(2:3)))
    v = NaN;
elseif isnan(vals(3))
    v = vals(2);
elseif isnan(vals(2))
    v = vals(3);
else
    v = mean(vals(2:3), 'omitnan');
end
end

function v = meanCnoPairRecruitedOnly(T, col)
if isempty(T) || ~ismember(col, T.Properties.VariableNames) || ~ismember('is_recruited_partition', T.Properties.VariableNames)
    v = NaN;
    return
end
mask = logical(double(T.is_recruited_partition));
if ~any(mask)
    v = NaN;
    return
end
T = T(mask, :);
v = meanCnoPair(T, col);
end

function v = meanCnoPairRecruitedCountByEvent(T, col)
if isempty(T) || ~ismember(col, T.Properties.VariableNames) || ~ismember('event_uid', T.Properties.VariableNames)
    v = NaN;
    return
end
vals = nan(1,3);
for e = 1:3
    Te = T(T.epoch_idx == e, {'event_uid', col});
    if isempty(Te)
        continue
    end
    [G,~] = findgroups(string(Te.event_uid));
    ev = splitapply(@(x) sum(logical(x(isfinite(x)))), double(Te.(col)), G);
    ev = double(ev(:));
    ev = ev(isfinite(ev));
    if ~isempty(ev)
        vals(e) = mean(ev, 'omitnan');
    end
end
if all(isnan(vals(2:3)))
    v = NaN;
elseif isnan(vals(3))
    v = vals(2);
elseif isnan(vals(2))
    v = vals(3);
else
    v = mean(vals(2:3), 'omitnan');
end
end


function ids = defaultSpiderCellMetricIDs()
ids = string({ ...
    'cell.frac_multicell_events', ...
    'cell.mean_cells_per_multicell_event', ...
    'cell.frac_soma_events', ...
    'cell.frac_toward_events', ...
    'cell.frac_away_events', ...
    'cell.frac_traverse_events', ...
    'cell.frac_net_toward_soma_events', ...
    'cell.frac_net_away_soma_events', ...
    'cell.mean_soma_net_distance_change_um', ...
    'cell.mean_soma_centroid_dist_um', ...
    'cell.mean_soma_start_dist_um', ...
    'cell.mean_soma_end_dist_um', ...
    'cell.mean_perimeter_um', ...
    'cell.mean_circularity', ...
    'cell.frac_propagative_events', ...
    'cell.mean_propagation_distance_um', ...
    'cell.mean_frac_outside_origin_cell', ...
    'cell.frac_escape_events_gt_threshold', ...
    'cell.frac_events_with_neighbor_recruitment', ...
    'cell.mean_neighbor_entry_delay_sec', ...
    'cell.mean_neighbor_peak_overlap_frac_of_cell', ...
    'cell.mean_neighbor_partitions_per_event'});
end

function ids = defaultSpiderSliceMetricIDs()
ids = string({ ...
    'slice.frac_active_cells', ...
    'slice.coincidence_rate', ...
    'slice.burst_fraction', ...
    'slice.mean_pairwise_cell_correlation'});
end

function cells = getEntityValueCells(polData, fieldName, nLabels)
if isfield(polData, fieldName) && iscell(polData.(fieldName)) && numel(polData.(fieldName)) == nLabels
    cells = polData.(fieldName)(:);
else
    cells = cell(nLabels, 1);
end
end

function key = chooseEntityKey(T, kind)
key = '';
candidates = {};
switch char(string(kind))
    case 'cell'
        candidates = {'primary_cell_uid', 'cell_uid', 'cell_id'};
    case 'event'
        candidates = {'event_uid', 'local_event_id'};
end
for i = 1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        key = candidates{i};
        return
    end
end
end

function vals = combinedPostEntityValues(T, valCol, entityKey)
vals = [];
if isempty(entityKey) || ~ismember(valCol, T.Properties.VariableNames) || ~ismember('epoch_idx', T.Properties.VariableNames)
    return
end
mask = T.epoch_idx >= 2 & isfinite(double(T.(valCol)));
if ~any(mask)
    return
end
Tp = T(mask, :);
entityID = makeEntityID(Tp, entityKey);
if isempty(entityID)
    return
end
G = findgroups(entityID);
vals = splitapply(@(x) mean(double(x), 'omitnan'), Tp.(valCol), G);
vals = double(vals(:));
vals = vals(isfinite(vals));
end

function vals = combinedPostValuesByEventAny(T, valCol)
vals = [];
if ~ismember('event_uid', T.Properties.VariableNames) || ~ismember(valCol, T.Properties.VariableNames) || ~ismember('epoch_idx', T.Properties.VariableNames)
    return
end
mask = T.epoch_idx >= 2;
if ~any(mask)
    return
end
Tp = T(mask, {'file_id','event_uid', valCol});
G = findgroups(makeEntityID(Tp, 'event_uid'));
vals = splitapply(@(x) any(logical(x(isfinite(x)))), double(Tp.(valCol)), G);
vals = double(vals(:));
vals = vals(isfinite(vals));
end

function vals = combinedPostValuesByEventMean(T, valCol, maskCol, maskValue)
vals = [];
if ~ismember('event_uid', T.Properties.VariableNames) || ~ismember(valCol, T.Properties.VariableNames) || ~ismember('epoch_idx', T.Properties.VariableNames)
    return
end
mask = T.epoch_idx >= 2;
if nargin >= 3 && ~isempty(maskCol) && ismember(maskCol, T.Properties.VariableNames)
    mask = mask & logical(double(T.(maskCol)) == double(maskValue));
end
if ~any(mask)
    return
end
Tp = T(mask, {'file_id','event_uid', valCol});
G = findgroups(makeEntityID(Tp, 'event_uid'));
vals = splitapply(@(x) mean(double(x), 'omitnan'), Tp.(valCol), G);
vals = double(vals(:));
vals = vals(isfinite(vals));
end

function vals = combinedPostValuesByEventCount(T, valCol)
vals = [];
if ~ismember('event_uid', T.Properties.VariableNames) || ~ismember(valCol, T.Properties.VariableNames) || ~ismember('epoch_idx', T.Properties.VariableNames)
    return
end
mask = T.epoch_idx >= 2;
if ~any(mask)
    return
end
Tp = T(mask, {'file_id','event_uid', valCol});
G = findgroups(makeEntityID(Tp, 'event_uid'));
vals = splitapply(@(x) sum(logical(x(isfinite(x)))), double(Tp.(valCol)), G);
vals = double(vals(:));
vals = vals(isfinite(vals));
end

function id = makeEntityID(T, key)
id = [];
if ~ismember(key, T.Properties.VariableNames)
    return
end
if ismember('file_id', T.Properties.VariableNames)
    id = string(T.file_id) + "|" + string(T.(key));
else
    id = string(T.(key));
end
end
