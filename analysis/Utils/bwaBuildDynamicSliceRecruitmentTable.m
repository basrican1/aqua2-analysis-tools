function T = bwaBuildDynamicSliceRecruitmentTable(result, varargin)
% Build promoted slice recruitment summaries on demand from event partitions.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'result', @isstruct);
addParameter(p, 'transform', 'delta', @(x) ischar(x) || isstring(x));
parse(p, result, varargin{:});

T = table();
if ~(isfield(result, 'metrics') && isstruct(result.metrics) && ...
        isfield(result.metrics, 'tables') && isstruct(result.metrics.tables) && ...
        isfield(result.metrics.tables, 'event_partition'))
    return
end

eventPartition = result.metrics.tables.event_partition;
if ~istable(eventPartition) || isempty(eventPartition)
    return
end

need = {'event_uid','animal_id','slice_num','file_id','epoch_idx','epoch_name'};
if ~all(ismember(need, eventPartition.Properties.VariableNames))
    return
end

P = eventPartition;
if ~ismember('is_recruited_partition', P.Properties.VariableNames)
    return
end
P.event_uid = string(P.event_uid);
P.animal_id = string(P.animal_id);
P.slice_num = string(P.slice_num);
P.file_id = string(P.file_id);
P.epoch_name = string(P.epoch_name);

if ~ismember('cell_entry_delay_from_primary_sec', P.Properties.VariableNames)
    if all(ismember({'cell_entry_delay_from_primary_frames','frame_rate_s'}, P.Properties.VariableNames))
        P.cell_entry_delay_from_primary_sec = double(P.cell_entry_delay_from_primary_frames) .* double(P.frame_rate_s);
    else
        P.cell_entry_delay_from_primary_sec = NaN(height(P),1);
    end
end
if ~ismember('cell_peak_overlap_frac_of_cell', P.Properties.VariableNames)
    P.cell_peak_overlap_frac_of_cell = NaN(height(P),1);
end

primaryEvents = unique(P(:, {'event_uid','animal_id','slice_num','file_id','epoch_idx','epoch_name'}), 'rows', 'stable');
if isempty(primaryEvents)
    return
end

nonPrimary = P(logical(P.is_recruited_partition), :);
primaryEvents.n_neighbor_partitions = zeros(height(primaryEvents),1);
primaryEvents.has_neighbor_recruitment = false(height(primaryEvents),1);
primaryEvents.mean_neighbor_entry_delay_sec = NaN(height(primaryEvents),1);
primaryEvents.mean_neighbor_peak_overlap_frac_of_cell = NaN(height(primaryEvents),1);

if ~isempty(nonPrimary)
    [G, grp] = localGroupRows(nonPrimary, {'event_uid'});
    nG = max([G;0]);
    perEvent = grp;
    perEvent.n_neighbor_partitions = zeros(nG,1);
    perEvent.mean_neighbor_entry_delay_sec = NaN(nG,1);
    perEvent.mean_neighbor_peak_overlap_frac_of_cell = NaN(nG,1);
    for k = 1:nG
        Tk = nonPrimary(G == k,:);
        perEvent.n_neighbor_partitions(k) = height(Tk);
        x = double(Tk.cell_entry_delay_from_primary_sec);
        x = x(isfinite(x));
        if ~isempty(x)
            perEvent.mean_neighbor_entry_delay_sec(k) = mean(x, 'omitnan');
        end
        x = double(Tk.cell_peak_overlap_frac_of_cell);
        x = x(isfinite(x));
        if ~isempty(x)
            perEvent.mean_neighbor_peak_overlap_frac_of_cell(k) = mean(x, 'omitnan');
        end
    end
    perEvent.has_neighbor_recruitment = perEvent.n_neighbor_partitions > 0;
    primaryEvents = outerjoin(primaryEvents, perEvent, 'Keys', {'event_uid'}, 'Type', 'left', 'MergeKeys', true);
    primaryEvents = coalesceJoinedColumns(primaryEvents);
end

primaryEvents.n_neighbor_partitions(~isfinite(primaryEvents.n_neighbor_partitions)) = 0;
primaryEvents.has_neighbor_recruitment = logical(primaryEvents.has_neighbor_recruitment);

[G, grp] = localGroupRows(primaryEvents, {'animal_id','slice_num','file_id','epoch_idx','epoch_name'});
nG = max([G;0]);
if nG == 0
    T = table();
    return
end

T = grp;
T.n_primary_events = zeros(nG,1);
T.n_primary_events_with_neighbor_recruitment = zeros(nG,1);
T.frac_primary_events_with_neighbor_recruitment = NaN(nG,1);
T.mean_neighbor_entry_delay_sec = NaN(nG,1);
T.mean_neighbor_peak_overlap_frac_of_cell = NaN(nG,1);
T.mean_neighbor_partitions_per_event = NaN(nG,1);

for k = 1:nG
    Tk = primaryEvents(G == k,:);
    T.n_primary_events(k) = height(Tk);
    T.n_primary_events_with_neighbor_recruitment(k) = nnz(Tk.has_neighbor_recruitment);
    if height(Tk) > 0
        T.frac_primary_events_with_neighbor_recruitment(k) = mean(double(Tk.has_neighbor_recruitment), 'omitnan');
        x = double(Tk.mean_neighbor_entry_delay_sec); x = x(isfinite(x));
        if ~isempty(x), T.mean_neighbor_entry_delay_sec(k) = mean(x, 'omitnan'); end
        x = double(Tk.mean_neighbor_peak_overlap_frac_of_cell); x = x(isfinite(x));
        if ~isempty(x), T.mean_neighbor_peak_overlap_frac_of_cell(k) = mean(x, 'omitnan'); end
        x = double(Tk.n_neighbor_partitions); x = x(isfinite(x));
        if ~isempty(x), T.mean_neighbor_partitions_per_event(k) = mean(x, 'omitnan'); end
    end
end

transform = lower(char(string(p.Results.transform)));
T.frac_primary_events_with_neighbor_recruitment_norm = localNormByBaseline(T, ...
    'frac_primary_events_with_neighbor_recruitment', {'animal_id','slice_num','file_id'}, 1e-6, transform);
T.mean_neighbor_entry_delay_sec_norm = localNormByBaseline(T, ...
    'mean_neighbor_entry_delay_sec', {'animal_id','slice_num','file_id'}, 1e-6, transform);
T.mean_neighbor_peak_overlap_frac_of_cell_norm = localNormByBaseline(T, ...
    'mean_neighbor_peak_overlap_frac_of_cell', {'animal_id','slice_num','file_id'}, 1e-6, transform);
T.mean_neighbor_partitions_per_event_norm = localNormByBaseline(T, ...
    'mean_neighbor_partitions_per_event', {'animal_id','slice_num','file_id'}, 1e-6, transform);
end

function T = coalesceJoinedColumns(T)
vars = string(T.Properties.VariableNames);
for i = 1:numel(vars)
    v = vars(i);
    if endsWith(v, "_primaryEvents")
        base = extractBefore(v, strlength(v) - strlength("_primaryEvents") + 1);
        alt = base + "_perEvent";
        if any(vars == alt)
            if isnumeric(T.(char(v))) || islogical(T.(char(v)))
                x = T.(char(v));
                y = T.(char(alt));
                miss = ~isfinite(double(x));
                x(miss) = y(miss);
                T.(char(base)) = x;
            else
                x = string(T.(char(v)));
                y = string(T.(char(alt)));
                miss = strlength(x) == 0;
                x(miss) = y(miss);
                T.(char(base)) = x;
            end
        end
    end
end
drop = vars(endsWith(vars, "_primaryEvents") | endsWith(vars, "_perEvent"));
drop = drop(ismember(drop, string(T.Properties.VariableNames)));
if ~isempty(drop)
    T = removevars(T, cellstr(drop));
end
end

function normVals = localNormByBaseline(T, metricCol, keyCols, epsVal, transform)
normVals = NaN(height(T), 1);
if ~ismember(metricCol, T.Properties.VariableNames) || ~ismember('epoch_idx', T.Properties.VariableNames)
    return
end
if nargin < 4 || isempty(epsVal)
    epsVal = 0;
end
if nargin < 5 || isempty(transform)
    transform = 'delta';
end
keyCols = keyCols(ismember(keyCols, T.Properties.VariableNames));
if isempty(keyCols)
    return
end
[G, ~] = localGroupRows(T, keyCols);
nG = max([G;0]);
for k = 1:nG
    idx = (G == k);
    Tk = T(idx, :);
    b = double(Tk.(metricCol)(Tk.epoch_idx == 1));
    b = b(isfinite(b));
    if isempty(b)
        continue
    end
    x = double(Tk.(metricCol));
    if strcmp(transform, 'log2_fold')
        x = x + epsVal;
        b0 = bwaComputeBaselineValue(x, Tk.epoch_idx == 1, transform);
    else
        b0 = bwaComputeBaselineValue(x, Tk.epoch_idx == 1, transform);
    end
    base = repmat(b0, size(x));
    normVals(idx) = bwaApplyTransform(x, base, transform);
end
end

function [G, grpTbl] = localGroupRows(T, cols)
n = height(T);
if n == 0
    G = zeros(0,1);
    grpTbl = T(:, cols);
    return
end
K = strings(n,1);
for i = 1:numel(cols)
    K = K + "|" + toKey(T.(cols{i}));
end
[~,~,G] = unique(K, 'stable');
G = full(double(G(:)));
firstIdx = accumarray(G, (1:n)', [], @min);
grpTbl = T(firstIdx, cols);
end

function s = toKey(v)
if isstring(v)
    s = v(:);
elseif iscell(v)
    s = strings(numel(v),1);
    for i = 1:numel(v)
        s(i) = string(v{i});
    end
elseif isnumeric(v) || islogical(v)
    s = string(v(:));
else
    s = string(v(:));
end
end
