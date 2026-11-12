function M = build_local_metrics_schema(canonicalOut, canonicalNorm, ctx, epochSpec, varargin)
% Build unified table-native local metrics schema.
%
% Output:
%   M.tables.event_canonical
%   M.tables.event_spread
%   M.tables.cell_spread
%   M.tables.cell_frequency
%   M.tables.cell_canonical
%   M.tables.event_soma
%   M.tables.cell_soma
%   M.tables.event_shape
%   M.tables.cell_shape
%   M.tables.event_propagation
%   M.tables.cell_propagation
%   M.tables.event_boundary
%   M.tables.cell_boundary
%   M.tables.event_partition
%   M.tables.cell_recruitment
%   M.tables.slice_network
%   M.registry

p = inputParser;
p.StructExpand = false;
addRequired(p, 'canonicalOut', @isstruct);
addRequired(p, 'canonicalNorm', @isstruct);
addRequired(p, 'ctx', @isstruct);
addRequired(p, 'epochSpec', @isstruct);
addParameter(p, 'spreadThresholdPct', 10, @isnumeric);
addParameter(p, 'propagationThresholdPx', 5, @isnumeric);
addParameter(p, 'escapeThresholdFrac', 0.10, @isnumeric);
addParameter(p, 'coincidenceWindowSec', 3, @isnumeric);
addParameter(p, 'burstBinSec', 5, @isnumeric);
addParameter(p, 'burstMinEvents', 3, @isnumeric);
addParameter(p, 'corrBinSec', 5, @isnumeric);
addParameter(p, 'somaNetMoveThresholdPx', 2, @isnumeric);
addParameter(p, 'normalizationTransform', 'delta', @(x) ischar(x) || isstring(x));
parse(p, canonicalOut, canonicalNorm, ctx, epochSpec, varargin{:});
localMetricNormTransform(p.Results.normalizationTransform);

M = struct();
M.tables = struct();
M.tables.event_canonical = table();
M.tables.event_spread = table();
M.tables.cell_spread = table();
M.tables.cell_frequency = table();
M.tables.cell_canonical = table();
M.tables.event_soma = table();
M.tables.cell_soma = table();
M.tables.event_shape = table();
M.tables.cell_shape = table();
M.tables.event_propagation = table();
M.tables.cell_propagation = table();
M.tables.event_boundary = table();
M.tables.cell_boundary = table();
M.tables.event_partition = table();
M.tables.cell_recruitment = table();
M.tables.slice_network = table();
M.registry = emptyRegistry();

% 1) Cell frequency table from canonical frequency output.
if isfield(canonicalOut, 'frequency') && isstruct(canonicalOut.frequency) && ...
        isfield(canonicalOut.frequency, 'ratesByCellEpoch') && istable(canonicalOut.frequency.ratesByCellEpoch)
    M.tables.cell_frequency = buildCellFrequencyTable(canonicalOut, canonicalNorm);
    M.registry = addRegistryRow(M.registry, ...
        "cell.frequency", "cell", "cell_frequency", ...
        "rate_events_per_min", "rate_events_per_min_norm", true, "by_cell", "summary_line");
end

% 2) Spread metrics from event overlap traces.
try
    [eventSpread, cellSpread] = computeCellSpreadTables(canonicalOut, ctx, epochSpec, p.Results.spreadThresholdPct);
    M.tables.event_spread = eventSpread;
    M.tables.cell_spread = cellSpread;
    if ~isempty(cellSpread)
        M.registry = addRegistryRow(M.registry, ...
            "cell.frac_multicell_events", "cell", "cell_spread", ...
            "frac_multicell_events", "frac_multicell_events_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_cells_per_multicell_event", "cell", "cell_spread", ...
            "mean_cells_per_multicell_event", "mean_cells_per_multicell_event_norm", true, "by_cell", "summary_line");
    end
catch ME
    warning('build_local_metrics_schema:SpreadBuildFailed', ...
        'Failed building spread metrics tables: %s', ME.message);
end

% 2b) Promoted cell summaries for event-native canonical metrics.
try
    [eventCanonical, cellCanonical] = computeCanonicalTables(canonicalOut, canonicalNorm, epochSpec);
    M.tables.event_canonical = eventCanonical;
    M.tables.cell_canonical = cellCanonical;
    if ~isempty(eventCanonical)
        M.registry = addRegistryRow(M.registry, ...
            "event.area_um2", "event", "event_canonical", ...
            "area_um2", "area_um2_norm", true, "by_event", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "event.duration_sec", "event", "event_canonical", ...
            "duration_sec", "duration_sec_norm", true, "by_event", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "event.rise19_sec", "event", "event_canonical", ...
            "rise19_sec", "rise19_sec_norm", true, "by_event", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "event.fall91_sec", "event", "event_canonical", ...
            "fall91_sec", "fall91_sec_norm", true, "by_event", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "event.decay_tau_sec", "event", "event_canonical", ...
            "decay_tau_sec", "decay_tau_sec_norm", true, "by_event", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "event.signal_auc_dff_sec", "event", "event_canonical", ...
            "dff_auc_sec", "dff_auc_sec_norm", true, "by_event", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "event.peak_amplitude_dff", "event", "event_canonical", ...
            "dff_max2", "dff_max2_norm", true, "by_event", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "event.max_propagation_speed_um_per_sec", "event", "event_canonical", ...
            "max_prop_speed_um_per_sec", "max_prop_speed_um_per_sec_norm", true, "by_event", "summary_line");
    end
    if ~isempty(cellCanonical)
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_event_area_um2", "cell", "cell_canonical", ...
            "mean_event_area_um2", "mean_event_area_um2_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_event_duration_sec", "cell", "cell_canonical", ...
            "mean_event_duration_sec", "mean_event_duration_sec_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_event_rise19_sec", "cell", "cell_canonical", ...
            "mean_event_rise19_sec", "mean_event_rise19_sec_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_event_fall91_sec", "cell", "cell_canonical", ...
            "mean_event_fall91_sec", "mean_event_fall91_sec_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_event_decay_tau_sec", "cell", "cell_canonical", ...
            "mean_event_decay_tau_sec", "mean_event_decay_tau_sec_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_event_signal_auc_dff_sec", "cell", "cell_canonical", ...
            "mean_event_signal_auc_dff_sec", "mean_event_signal_auc_dff_sec_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_event_peak_amplitude_dff", "cell", "cell_canonical", ...
            "mean_event_peak_amplitude_dff", "mean_event_peak_amplitude_dff_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_event_max_propagation_speed_um_per_sec", "cell", "cell_canonical", ...
            "mean_event_max_propagation_speed_um_per_sec", "mean_event_max_propagation_speed_um_per_sec_norm", true, "by_cell", "summary_line");
    end
catch ME
    warning('build_local_metrics_schema:CanonicalMetricBuildFailed', ...
        'Failed building canonical event/cell metrics tables: %s', ME.message);
end

% 3) Soma/landmark metrics from canonical event_cell_map columns.
try
    [eventSoma, cellSoma] = computeSomaTables(canonicalOut, epochSpec, p.Results.somaNetMoveThresholdPx);
    M.tables.event_soma = eventSoma;
    M.tables.cell_soma = cellSoma;
    if ~isempty(cellSoma)
        M.registry = addRegistryRow(M.registry, ...
            "cell.frac_soma_events", "cell", "cell_soma", ...
            "frac_soma_events", "frac_soma_events_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.frac_toward_events", "cell", "cell_soma", ...
            "frac_toward_events", "frac_toward_events_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.frac_away_events", "cell", "cell_soma", ...
            "frac_away_events", "frac_away_events_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.frac_traverse_events", "cell", "cell_soma", ...
            "frac_traverse_events", "frac_traverse_events_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_motion_index", "cell", "cell_soma", ...
            "mean_motion_index", "", false, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_soma_centroid_dist_px", "cell", "cell_soma", ...
            "mean_soma_centroid_dist_px", "mean_soma_centroid_dist_px_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_soma_start_dist_px", "cell", "cell_soma", ...
            "mean_soma_start_dist_px", "mean_soma_start_dist_px_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_soma_end_dist_px", "cell", "cell_soma", ...
            "mean_soma_end_dist_px", "mean_soma_end_dist_px_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_soma_centroid_dist_um", "cell", "cell_soma", ...
            "mean_soma_centroid_dist_um", "mean_soma_centroid_dist_um_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_soma_start_dist_um", "cell", "cell_soma", ...
            "mean_soma_start_dist_um", "mean_soma_start_dist_um_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_soma_end_dist_um", "cell", "cell_soma", ...
            "mean_soma_end_dist_um", "mean_soma_end_dist_um_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.frac_net_toward_soma_events", "cell", "cell_soma", ...
            "frac_net_toward_soma_events", "frac_net_toward_soma_events_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.frac_net_away_soma_events", "cell", "cell_soma", ...
            "frac_net_away_soma_events", "frac_net_away_soma_events_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.frac_net_stationary_soma_events", "cell", "cell_soma", ...
            "frac_net_stationary_soma_events", "frac_net_stationary_soma_events_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_soma_net_distance_change_um", "cell", "cell_soma", ...
            "mean_soma_net_distance_change_um", "mean_soma_net_distance_change_um_norm", true, "by_cell", "summary_line");
    end
catch ME
    warning('build_local_metrics_schema:SomaBuildFailed', ...
        'Failed building soma metrics tables: %s', ME.message);
end

% 4) Shape metrics from canonical events columns.
try
    [eventShape, cellShape] = computeShapeTables(canonicalOut, canonicalNorm, epochSpec);
    M.tables.event_shape = eventShape;
    M.tables.cell_shape = cellShape;
    if ~isempty(cellShape)
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_perimeter_px", "cell", "cell_shape", ...
            "mean_perimeter_px", "mean_perimeter_px_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_circularity", "cell", "cell_shape", ...
            "mean_circularity", "mean_circularity_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_perimeter_um", "cell", "cell_shape", ...
            "mean_perimeter_um", "mean_perimeter_um_norm", true, "by_cell", "summary_line");
    end
catch ME
    warning('build_local_metrics_schema:ShapeBuildFailed', ...
        'Failed building shape metrics tables: %s', ME.message);
end

% 5) Propagation metrics from canonical events columns.
try
    [eventPropagation, cellPropagation] = computePropagationTables(canonicalOut, canonicalNorm, epochSpec, p.Results.propagationThresholdPx);
    M.tables.event_propagation = eventPropagation;
    M.tables.cell_propagation = cellPropagation;
    if ~isempty(cellPropagation)
        M.registry = addRegistryRow(M.registry, ...
            "cell.frac_propagative_events", "cell", "cell_propagation", ...
            "frac_propagative_events", "frac_propagative_events_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_propagation_distance_px", "cell", "cell_propagation", ...
            "mean_propagation_distance_px", "mean_propagation_distance_px_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_propagation_distance_um", "cell", "cell_propagation", ...
            "mean_propagation_distance_um", "mean_propagation_distance_um_norm", true, "by_cell", "summary_line");
    end
catch ME
    warning('build_local_metrics_schema:PropagationBuildFailed', ...
        'Failed building propagation metrics tables: %s', ME.message);
end

% 6) Boundary escape metrics from canonical event_cell_map primary rows.
try
    [eventBoundary, cellBoundary] = computeBoundaryTables(canonicalOut, canonicalNorm, epochSpec, p.Results.escapeThresholdFrac);
    M.tables.event_boundary = eventBoundary;
    M.tables.cell_boundary = cellBoundary;
    if ~isempty(cellBoundary)
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_frac_outside_origin_cell", "cell", "cell_boundary", ...
            "mean_frac_outside_origin_cell", "mean_frac_outside_origin_cell_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.frac_escape_events_gt_threshold", "cell", "cell_boundary", ...
            "frac_escape_events_gt_threshold", "frac_escape_events_gt_threshold_norm", true, "by_cell", "summary_line");
    end
catch ME
    warning('build_local_metrics_schema:BoundaryBuildFailed', ...
        'Failed building boundary metrics tables: %s', ME.message);
end

% 7) Event-cell partition / recruitment metrics from canonical event_cell_map rows.
try
    [eventPartition, cellRecruitment, ~] = computePartitionTables(canonicalOut, epochSpec);
    M.tables.event_partition = eventPartition;
    M.tables.cell_recruitment = cellRecruitment;
    if ~isempty(cellRecruitment)
        M.registry = addRegistryRow(M.registry, ...
            "cell.frac_events_with_neighbor_recruitment", "cell", "cell_recruitment", ...
            "frac_events_with_neighbor_recruitment", "frac_events_with_neighbor_recruitment_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_neighbor_entry_delay_frames", "cell", "cell_recruitment", ...
            "mean_neighbor_entry_delay_frames", "mean_neighbor_entry_delay_frames_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_neighbor_peak_overlap_frac_of_cell", "cell", "cell_recruitment", ...
            "mean_neighbor_peak_overlap_frac_of_cell", "mean_neighbor_peak_overlap_frac_of_cell_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_neighbor_partitions_per_event", "cell", "cell_recruitment", ...
            "mean_neighbor_partitions_per_event", "mean_neighbor_partitions_per_event_norm", true, "by_cell", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "cell.mean_neighbor_entry_delay_sec", "cell", "cell_recruitment", ...
            "mean_neighbor_entry_delay_sec", "mean_neighbor_entry_delay_sec_norm", true, "by_cell", "summary_line");
    end
catch ME
    warning('build_local_metrics_schema:PartitionBuildFailed', ...
        'Failed building partition metrics tables: %s', ME.message);
end

% 8) Slice-level network/synchrony metrics from primary-assigned events.
try
    sliceNetwork = computeSliceNetworkTables(canonicalOut, canonicalNorm, epochSpec, ...
        p.Results.coincidenceWindowSec, p.Results.burstBinSec, p.Results.burstMinEvents, p.Results.corrBinSec);
    M.tables.slice_network = sliceNetwork;
    if ~isempty(sliceNetwork)
        M.registry = addRegistryRow(M.registry, ...
            "slice.frac_active_cells", "slice", "slice_network", ...
            "frac_active_cells", "frac_active_cells_norm", true, "by_slice", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "slice.coincidence_rate", "slice", "slice_network", ...
            "coincidence_rate", "coincidence_rate_norm", true, "by_slice", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "slice.burst_fraction", "slice", "slice_network", ...
            "burst_fraction", "burst_fraction_norm", true, "by_slice", "summary_line");
        M.registry = addRegistryRow(M.registry, ...
            "slice.mean_pairwise_cell_correlation", "slice", "slice_network", ...
            "mean_pairwise_cell_correlation", "mean_pairwise_cell_correlation_norm", true, "by_slice", "summary_line");
    end
catch ME
    warning('build_local_metrics_schema:SliceNetworkBuildFailed', ...
        'Failed building slice network metrics tables: %s', ME.message);
end

M.registry = bwaEnrichMetricRegistry(M.registry);
end

function T = buildCellFrequencyTable(canonicalOut, canonicalNorm)
T = canonicalOut.frequency.ratesByCellEpoch;
if ismember('cell_uid', T.Properties.VariableNames) && ~ismember('primary_cell_uid', T.Properties.VariableNames)
    T.primary_cell_uid = string(T.cell_uid);
end
if ~ismember('primary_cell_uid', T.Properties.VariableNames)
    T.primary_cell_uid = repmat("", height(T), 1);
end
T.primary_cell_uid = string(T.primary_cell_uid);

T.rate_events_per_min_norm = NaN(height(T), 1);
if isfield(canonicalNorm, 'frequency') && isstruct(canonicalNorm.frequency) && ...
        isfield(canonicalNorm.frequency, 'ratesByCellEpoch') && istable(canonicalNorm.frequency.ratesByCellEpoch)
    TN = canonicalNorm.frequency.ratesByCellEpoch;
    if ismember('cell_uid', TN.Properties.VariableNames) && ~ismember('primary_cell_uid', TN.Properties.VariableNames)
        TN.primary_cell_uid = string(TN.cell_uid);
    end
    if all(ismember({'file_id','epoch_idx','primary_cell_uid','rate_events_per_min'}, TN.Properties.VariableNames))
        K = {'file_id','epoch_idx','primary_cell_uid'};
        TL = outerjoin(T(:, K), TN(:, [K {'rate_events_per_min'}]), ...
            'Keys', K, 'Type', 'left', 'MergeKeys', true);
        if ismember('rate_events_per_min_TN', TL.Properties.VariableNames)
            T.rate_events_per_min_norm = TL.rate_events_per_min_TN;
        elseif ismember('rate_events_per_min', TL.Properties.VariableNames)
            T.rate_events_per_min_norm = TL.rate_events_per_min;
        end
    end
end
end

function [eventSpread, cellSpread] = computeCellSpreadTables(canonicalOut, ctx, epochSpec, thresholdPct)
eventSpread = table();
cellSpread = table();
if ~(isfield(canonicalOut, 'subset') && isstruct(canonicalOut.subset) && ...
        isfield(canonicalOut.subset, 'events') && istable(canonicalOut.subset.events))
    return
end
if ~isfield(ctx, 'associatedCoverageTraces') || isempty(ctx.associatedCoverageTraces)
    return
end

E = canonicalOut.subset.events;
if ~ismember('epoch_idx', E.Properties.VariableNames) || ~ismember('epoch_name', E.Properties.VariableNames)
    A = assignEpochsToEvents(E, epochSpec, 'method', 'auto');
    E = A.events;
end
if ~ismember('primary_cell_uid', E.Properties.VariableNames)
    E.primary_cell_uid = repmat("", height(E), 1);
end
E = E(E.epoch_idx > 0, :);

legacySpread = legacySpreadLookup(ctx, thresholdPct);
if isempty(legacySpread)
    return
end
if ~all(ismember({'file_id','local_event_id'}, E.Properties.VariableNames))
    return
end

J = outerjoin(E, legacySpread, 'Keys', {'file_id','local_event_id'}, 'Type', 'left', 'MergeKeys', true);
if ~ismember('spread_count', J.Properties.VariableNames)
    J.spread_count = NaN(height(J),1);
end
J.spread_count = double(J.spread_count);
J.is_multicell = J.spread_count > 1 & isfinite(J.spread_count);

keepCols = {'event_uid','animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name','spread_count','is_multicell','spread_threshold_pct'};
keepCols = keepCols(ismember(keepCols, J.Properties.VariableNames));
eventSpread = J(:, keepCols);
if isempty(eventSpread)
    return
end

eventSpread.primary_cell_uid = string(eventSpread.primary_cell_uid);
eventSpread = eventSpread(string(eventSpread.primary_cell_uid) ~= "", :);
if isempty(eventSpread)
    return
end

eventSpread.spread_count_norm = localNormByCellBaseline(eventSpread, 'spread_count', ...
    {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);

keyCols = {'animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name'};
keyCols = keyCols(ismember(keyCols, eventSpread.Properties.VariableNames));
[G, grp] = localGroupRows(eventSpread, keyCols);
nG = max([G;0]);
if nG == 0
    return
end

nTotal = zeros(nG,1);
nMulti = zeros(nG,1);
fracMulti = NaN(nG,1);
meanSpread = NaN(nG,1);
nMissingSpread = zeros(nG,1);
for k = 1:nG
    idx = (G == k);
    x = double(eventSpread.spread_count(idx));
    good = isfinite(x);
    nTotal(k) = sum(good);
    nMissingSpread(k) = sum(~good);
    if nTotal(k) == 0
        continue
    end
    xm = x(good);
    nMulti(k) = sum(xm > 1);
    fracMulti(k) = nMulti(k) / nTotal(k);
    if nMulti(k) > 0
        meanSpread(k) = mean(xm(xm > 1), 'omitnan');
    end
end

cellSpread = grp;
cellSpread.n_total_events = nTotal;
cellSpread.n_multicell_events = nMulti;
cellSpread.n_missing_spread = nMissingSpread;
cellSpread.frac_multicell_events = fracMulti;
cellSpread.mean_cells_per_multicell_event = meanSpread;
cellSpread.spread_threshold_pct = repmat(double(thresholdPct), height(cellSpread), 1);
cellSpread.frac_multicell_events_norm = localNormByCellBaseline(cellSpread, ...
    'frac_multicell_events', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellSpread.mean_cells_per_multicell_event_norm = localNormByCellBaseline(cellSpread, ...
    'mean_cells_per_multicell_event', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
end

function [eventSoma, cellSoma] = computeSomaTables(canonicalOut, epochSpec, somaNetMoveThresholdPx)
eventSoma = table();
cellSoma = table();

if ~(isfield(canonicalOut, 'subset') && isstruct(canonicalOut.subset) && ...
        isfield(canonicalOut.subset, 'events') && istable(canonicalOut.subset.events) && ...
        isfield(canonicalOut.subset, 'event_cell_map') && istable(canonicalOut.subset.event_cell_map))
    return
end

E = canonicalOut.subset.events;
M = canonicalOut.subset.event_cell_map;

needMap = {'event_uid','cell_uid','is_primary_cell','soma_covered_ge50','motion_class','motion_index'};
optMap = {'soma_centroid_dist_px','soma_start_dist_px','soma_end_dist_px'};
if ~all(ismember(needMap, M.Properties.VariableNames))
    return
end

if ~ismember('epoch_idx', E.Properties.VariableNames) || ~ismember('epoch_name', E.Properties.VariableNames)
    A = assignEpochsToEvents(E, epochSpec, 'method', 'auto');
    E = A.events;
end

joinCols = {'event_uid','animal_id','slice_num','file_id','epoch_idx','epoch_name'};
joinCols = joinCols(ismember(joinCols, E.Properties.VariableNames));
if ~ismember('event_uid', joinCols)
    return
end
J = outerjoin(M, E(:, joinCols), 'Keys', {'event_uid'}, 'Type', 'left', 'MergeKeys', true);
if ~all(ismember({'file_id','epoch_idx','epoch_name'}, J.Properties.VariableNames))
    return
end

J = J(J.epoch_idx > 0, :);
if isempty(J)
    return
end

J.primary_cell_uid = string(J.cell_uid);
J.soma_covered_ge50 = logical(J.soma_covered_ge50);
J.motion_class = string(J.motion_class);
J.motion_index = double(J.motion_index);
for vn = string(optMap)
    if ismember(vn, string(J.Properties.VariableNames))
        J.(char(vn)) = double(J.(char(vn)));
    else
        J.(char(vn)) = NaN(height(J), 1);
    end
end

if ~ismember('animal_id', J.Properties.VariableNames)
    J.animal_id = repmat("", height(J), 1);
end
if ~ismember('slice_num', J.Properties.VariableNames)
    J.slice_num = repmat("", height(J), 1);
end

eventSoma = J(:, {'event_uid','animal_id','slice_num','file_id','primary_cell_uid', ...
    'epoch_idx','epoch_name','is_primary_cell','soma_covered_ge50','motion_class','motion_index', ...
    'soma_centroid_dist_px','soma_start_dist_px','soma_end_dist_px'});
eventSoma.soma_net_distance_change_px = eventSoma.soma_end_dist_px - eventSoma.soma_start_dist_px;
eventSoma.is_net_toward_soma = isfinite(eventSoma.soma_net_distance_change_px) & eventSoma.soma_net_distance_change_px <= -double(somaNetMoveThresholdPx);
eventSoma.is_net_away_soma = isfinite(eventSoma.soma_net_distance_change_px) & eventSoma.soma_net_distance_change_px >= double(somaNetMoveThresholdPx);
eventSoma.is_net_stationary_soma = isfinite(eventSoma.soma_net_distance_change_px) & abs(eventSoma.soma_net_distance_change_px) < double(somaNetMoveThresholdPx);
eventSoma = attachFileCalibration(eventSoma, canonicalOut);
if ismember('spatial_res', eventSoma.Properties.VariableNames)
    eventSoma.soma_centroid_dist_um = convertPxToUmColumn(eventSoma, 'soma_centroid_dist_px', 'spatial_res');
    eventSoma.soma_start_dist_um = convertPxToUmColumn(eventSoma, 'soma_start_dist_px', 'spatial_res');
    eventSoma.soma_end_dist_um = convertPxToUmColumn(eventSoma, 'soma_end_dist_px', 'spatial_res');
    eventSoma.soma_net_distance_change_um = convertPxToUmColumn(eventSoma, 'soma_net_distance_change_px', 'spatial_res');
    eventSoma.soma_centroid_dist_um_norm = localNormByCellBaseline(eventSoma, 'soma_centroid_dist_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
    eventSoma.soma_start_dist_um_norm = localNormByCellBaseline(eventSoma, 'soma_start_dist_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
    eventSoma.soma_end_dist_um_norm = localNormByCellBaseline(eventSoma, 'soma_end_dist_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
    eventSoma.soma_net_distance_change_um_norm = localNormByCellBaseline(eventSoma, 'soma_net_distance_change_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
end

keyCols = {'animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name'};
eventSomaPrimary = eventSoma(logical(eventSoma.is_primary_cell), :);
[G, grp] = localGroupRows(eventSoma, keyCols);
nG = max([G;0]);
if nG == 0
    return
end

nTotal = zeros(nG,1);
nSoma = zeros(nG,1);
nToward = zeros(nG,1);
nAway = zeros(nG,1);
nTraverse = zeros(nG,1);
nMotionValid = zeros(nG,1);
meanMotion = NaN(nG,1);
meanSomaCentroid = NaN(nG,1);
meanSomaStart = NaN(nG,1);
meanSomaEnd = NaN(nG,1);
fracNetToward = NaN(nG,1);
fracNetAway = NaN(nG,1);
fracNetStationary = NaN(nG,1);
meanNetChangePx = NaN(nG,1);

for k = 1:nG
    idx = (G == k);
    Tk = eventSoma(idx,:);
    nTotal(k) = height(Tk);
    if nTotal(k) == 0
        continue
    end
    nSoma(k) = nnz(Tk.soma_covered_ge50);

    mc = lower(string(Tk.motion_class));
    nToward(k) = nnz(mc == "toward");
    nAway(k) = nnz(mc == "away");
    nTraverse(k) = nnz(mc == "traverse");
    nMotionValid(k) = nToward(k) + nAway(k) + nTraverse(k);

    mi = double(Tk.motion_index);
    mi = mi(isfinite(mi));
    if ~isempty(mi)
        meanMotion(k) = mean(mi, 'omitnan');
    end

    ds = double(Tk.soma_centroid_dist_px);
    ds = ds(isfinite(ds));
    if ~isempty(ds)
        meanSomaCentroid(k) = mean(ds, 'omitnan');
    end
    ds = double(Tk.soma_start_dist_px);
    ds = ds(isfinite(ds));
    if ~isempty(ds)
        meanSomaStart(k) = mean(ds, 'omitnan');
    end
    ds = double(Tk.soma_end_dist_px);
    ds = ds(isfinite(ds));
    if ~isempty(ds)
        meanSomaEnd(k) = mean(ds, 'omitnan');
    end

    Tkp = eventSomaPrimary(string(eventSomaPrimary.animal_id) == string(grp.animal_id(k)) & ...
        string(eventSomaPrimary.slice_num) == string(grp.slice_num(k)) & ...
        string(eventSomaPrimary.file_id) == string(grp.file_id(k)) & ...
        string(eventSomaPrimary.primary_cell_uid) == string(grp.primary_cell_uid(k)) & ...
        double(eventSomaPrimary.epoch_idx) == double(grp.epoch_idx(k)), :);
    if ~isempty(Tkp)
        fracNetToward(k) = mean(double(Tkp.is_net_toward_soma), 'omitnan');
        fracNetAway(k) = mean(double(Tkp.is_net_away_soma), 'omitnan');
        fracNetStationary(k) = mean(double(Tkp.is_net_stationary_soma), 'omitnan');
        x = double(Tkp.soma_net_distance_change_px);
        x = x(isfinite(x));
        if ~isempty(x)
            meanNetChangePx(k) = mean(x, 'omitnan');
        end
    end
end

cellSoma = grp;
cellSoma.n_total_events = nTotal;
cellSoma.n_soma_events = nSoma;
cellSoma.frac_soma_events = nSoma ./ max(nTotal, 1);
cellSoma.n_toward_events = nToward;
cellSoma.n_away_events = nAway;
cellSoma.n_traverse_events = nTraverse;
cellSoma.n_motion_directional_events = nMotionValid;
cellSoma.frac_toward_events = nToward ./ max(nTotal, 1);
cellSoma.frac_away_events = nAway ./ max(nTotal, 1);
cellSoma.frac_traverse_events = nTraverse ./ max(nTotal, 1);
cellSoma.mean_motion_index = meanMotion;
cellSoma.mean_soma_centroid_dist_px = meanSomaCentroid;
cellSoma.mean_soma_start_dist_px = meanSomaStart;
cellSoma.mean_soma_end_dist_px = meanSomaEnd;
cellSoma.frac_net_toward_soma_events = fracNetToward;
cellSoma.frac_net_away_soma_events = fracNetAway;
cellSoma.frac_net_stationary_soma_events = fracNetStationary;
cellSoma.mean_soma_net_distance_change_px = meanNetChangePx;
cellSoma.soma_net_move_threshold_px = repmat(double(somaNetMoveThresholdPx), height(cellSoma), 1);

cellSoma.frac_soma_events_norm = localNormByCellBaseline(cellSoma, ...
    'frac_soma_events', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellSoma.frac_toward_events_norm = localNormByCellBaseline(cellSoma, ...
    'frac_toward_events', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellSoma.frac_away_events_norm = localNormByCellBaseline(cellSoma, ...
    'frac_away_events', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellSoma.frac_traverse_events_norm = localNormByCellBaseline(cellSoma, ...
    'frac_traverse_events', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellSoma.frac_net_toward_soma_events_norm = localNormByCellBaseline(cellSoma, ...
    'frac_net_toward_soma_events', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellSoma.frac_net_away_soma_events_norm = localNormByCellBaseline(cellSoma, ...
    'frac_net_away_soma_events', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellSoma.frac_net_stationary_soma_events_norm = localNormByCellBaseline(cellSoma, ...
    'frac_net_stationary_soma_events', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellSoma.mean_soma_centroid_dist_px_norm = localNormByCellBaseline(cellSoma, ...
    'mean_soma_centroid_dist_px', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
cellSoma.mean_soma_start_dist_px_norm = localNormByCellBaseline(cellSoma, ...
    'mean_soma_start_dist_px', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
cellSoma.mean_soma_end_dist_px_norm = localNormByCellBaseline(cellSoma, ...
    'mean_soma_end_dist_px', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
cellSoma.mean_soma_net_distance_change_px_norm = localNormByCellBaseline(cellSoma, ...
    'mean_soma_net_distance_change_px', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
cellSoma = attachFileCalibration(cellSoma, canonicalOut);
if ismember('spatial_res', cellSoma.Properties.VariableNames)
    cellSoma.mean_soma_centroid_dist_um = convertPxToUmColumn(cellSoma, 'mean_soma_centroid_dist_px', 'spatial_res');
    cellSoma.mean_soma_start_dist_um = convertPxToUmColumn(cellSoma, 'mean_soma_start_dist_px', 'spatial_res');
    cellSoma.mean_soma_end_dist_um = convertPxToUmColumn(cellSoma, 'mean_soma_end_dist_px', 'spatial_res');
    cellSoma.mean_soma_net_distance_change_um = convertPxToUmColumn(cellSoma, 'mean_soma_net_distance_change_px', 'spatial_res');
else
    cellSoma.mean_soma_centroid_dist_um = NaN(height(cellSoma),1);
    cellSoma.mean_soma_start_dist_um = NaN(height(cellSoma),1);
    cellSoma.mean_soma_end_dist_um = NaN(height(cellSoma),1);
    cellSoma.mean_soma_net_distance_change_um = NaN(height(cellSoma),1);
end
cellSoma.mean_soma_centroid_dist_um_norm = localNormByCellBaseline(cellSoma, ...
    'mean_soma_centroid_dist_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
cellSoma.mean_soma_start_dist_um_norm = localNormByCellBaseline(cellSoma, ...
    'mean_soma_start_dist_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
cellSoma.mean_soma_end_dist_um_norm = localNormByCellBaseline(cellSoma, ...
    'mean_soma_end_dist_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
cellSoma.mean_soma_net_distance_change_um_norm = localNormByCellBaseline(cellSoma, ...
    'mean_soma_net_distance_change_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
end

function [eventShape, cellShape] = computeShapeTables(canonicalOut, canonicalNorm, epochSpec)
eventShape = table();
cellShape = table();
need = {'event_uid','animal_id','slice_num','file_id','perimeter_px','circularity'};
[Eraw, Enorm] = canonicalEventTables(canonicalOut, canonicalNorm, epochSpec, need);
if isempty(Eraw)
    return
end

eventShape = Eraw(:, intersect({'event_uid','animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name','perimeter_px','circularity'}, Eraw.Properties.VariableNames));
if ~isempty(Enorm)
    eventShape = addNormColumns(eventShape, Enorm, {'perimeter_px','circularity'});
else
    eventShape.perimeter_px_norm = localNormByCellBaseline(eventShape, 'perimeter_px', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
    eventShape.circularity_norm = localNormByCellBaseline(eventShape, 'circularity', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
end
eventShape = attachFileCalibration(eventShape, canonicalOut);
if ismember('spatial_res', eventShape.Properties.VariableNames)
    eventShape.perimeter_um = convertPxToUmColumn(eventShape, 'perimeter_px', 'spatial_res');
    eventShape.perimeter_um_norm = localNormByCellBaseline(eventShape, 'perimeter_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
end

cellShape = summarizeCellMetricTable(eventShape, ...
    {'perimeter_px','mean_perimeter_px'; 'circularity','mean_circularity'});
cellShape.mean_perimeter_px_norm = localNormByCellBaseline(cellShape, 'mean_perimeter_px', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
cellShape.mean_circularity_norm = localNormByCellBaseline(cellShape, 'mean_circularity', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
cellShape = attachFileCalibration(cellShape, canonicalOut);
if ismember('spatial_res', cellShape.Properties.VariableNames)
    cellShape.mean_perimeter_um = convertPxToUmColumn(cellShape, 'mean_perimeter_px', 'spatial_res');
    cellShape.mean_perimeter_um_norm = localNormByCellBaseline(cellShape, 'mean_perimeter_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
end
end

function [eventPropagation, cellPropagation] = computePropagationTables(canonicalOut, canonicalNorm, epochSpec, thresholdPx)
eventPropagation = table();
cellPropagation = table();
need = {'event_uid','animal_id','slice_num','file_id','propagation_distance_px','propagation_angle_deg'};
[Eraw, Enorm] = canonicalEventTables(canonicalOut, canonicalNorm, epochSpec, need);
if isempty(Eraw)
    return
end

eventPropagation = Eraw(:, intersect({'event_uid','animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name','propagation_distance_px','propagation_angle_deg'}, Eraw.Properties.VariableNames));
eventPropagation.is_propagative = isfinite(eventPropagation.propagation_distance_px) & eventPropagation.propagation_distance_px >= double(thresholdPx);
if ~isempty(Enorm)
    eventPropagation = addNormColumns(eventPropagation, Enorm, {'propagation_distance_px'});
else
    eventPropagation.propagation_distance_px_norm = localNormByCellBaseline(eventPropagation, 'propagation_distance_px', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
end
eventPropagation.propagation_threshold_px = repmat(double(thresholdPx), height(eventPropagation), 1);
eventPropagation = attachFileCalibration(eventPropagation, canonicalOut);
if ismember('spatial_res', eventPropagation.Properties.VariableNames)
    eventPropagation.propagation_distance_um = convertPxToUmColumn(eventPropagation, 'propagation_distance_px', 'spatial_res');
    eventPropagation.propagation_distance_um_norm = localNormByCellBaseline(eventPropagation, 'propagation_distance_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
end

cellPropagation = summarizeCellMetricTable(eventPropagation, ...
    {'propagation_distance_px','mean_propagation_distance_px'});
if isempty(cellPropagation)
    return
end
[Gcnt, ~] = localGroupRows(eventPropagation, {'animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name'});
cellPropagation.n_total_events = accumarray(Gcnt, 1, [], @sum, 0);
cellPropagation.frac_propagative_events = localGroupedFraction(eventPropagation, {'animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name'}, 'is_propagative');
cellPropagation.propagation_threshold_px = repmat(double(thresholdPx), height(cellPropagation), 1);
cellPropagation.mean_propagation_distance_px_norm = localNormByCellBaseline(cellPropagation, 'mean_propagation_distance_px', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
cellPropagation.frac_propagative_events_norm = localNormByCellBaseline(cellPropagation, 'frac_propagative_events', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellPropagation = attachFileCalibration(cellPropagation, canonicalOut);
if ismember('spatial_res', cellPropagation.Properties.VariableNames)
    cellPropagation.mean_propagation_distance_um = convertPxToUmColumn(cellPropagation, 'mean_propagation_distance_px', 'spatial_res');
    cellPropagation.mean_propagation_distance_um_norm = localNormByCellBaseline(cellPropagation, 'mean_propagation_distance_um', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
end
end

function [eventCanonical, cellCanonical] = computeCanonicalTables(canonicalOut, canonicalNorm, epochSpec)
eventCanonical = table();
cellCanonical = table();
need = {'event_uid','animal_id','slice_num','file_id','area_um2','duration_sec','rise19_sec','fall91_sec', ...
    'decay_tau_sec','dff_auc_sec','dff_max2','max_prop_speed_um_per_sec'};
[Eraw, Enorm] = canonicalEventTables(canonicalOut, canonicalNorm, epochSpec, need);
if isempty(Eraw)
    return
end

eventCanonical = Eraw(:, intersect({'event_uid','animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name', ...
        'area_px2','area_um2','duration_sec','rise19_sec','fall91_sec','decay_tau_sec','dff_auc_sec','dff_max2','max_prop_speed_um_per_sec'}, ...
    Eraw.Properties.VariableNames));

if ~isempty(Enorm)
    eventCanonical = addNormColumns(eventCanonical, Enorm, {'area_um2','duration_sec','rise19_sec','fall91_sec', ...
        'decay_tau_sec','dff_auc_sec','dff_max2','max_prop_speed_um_per_sec'});
else
    keyCols = {'animal_id','slice_num','file_id','primary_cell_uid'};
    eventCanonical.area_um2_norm = localNormByCellBaseline(eventCanonical, 'area_um2', keyCols, 0);
    eventCanonical.duration_sec_norm = localNormByCellBaseline(eventCanonical, 'duration_sec', keyCols, 0);
    eventCanonical.rise19_sec_norm = localNormByCellBaseline(eventCanonical, 'rise19_sec', keyCols, 0);
    eventCanonical.fall91_sec_norm = localNormByCellBaseline(eventCanonical, 'fall91_sec', keyCols, 0);
    eventCanonical.decay_tau_sec_norm = localNormByCellBaseline(eventCanonical, 'decay_tau_sec', keyCols, 0);
    eventCanonical.dff_auc_sec_norm = localNormByCellBaseline(eventCanonical, 'dff_auc_sec', keyCols, 0);
    eventCanonical.dff_max2_norm = localNormByCellBaseline(eventCanonical, 'dff_max2', keyCols, 0);
    eventCanonical.max_prop_speed_um_per_sec_norm = localNormByCellBaseline(eventCanonical, 'max_prop_speed_um_per_sec', keyCols, 0);
end

cellCanonical = summarizeCellMetricTable(eventCanonical, { ...
    'area_um2','mean_event_area_um2'; ...
    'duration_sec','mean_event_duration_sec'; ...
    'rise19_sec','mean_event_rise19_sec'; ...
    'fall91_sec','mean_event_fall91_sec'; ...
    'decay_tau_sec','mean_event_decay_tau_sec'; ...
    'dff_auc_sec','mean_event_signal_auc_dff_sec'; ...
    'dff_max2','mean_event_peak_amplitude_dff'; ...
    'max_prop_speed_um_per_sec','mean_event_max_propagation_speed_um_per_sec'});

if isempty(cellCanonical)
    return
end

keyCols = {'animal_id','slice_num','file_id','primary_cell_uid'};
cellCanonical.mean_event_area_um2_norm = localNormByCellBaseline(cellCanonical, 'mean_event_area_um2', keyCols, 0);
cellCanonical.mean_event_duration_sec_norm = localNormByCellBaseline(cellCanonical, 'mean_event_duration_sec', keyCols, 0);
cellCanonical.mean_event_rise19_sec_norm = localNormByCellBaseline(cellCanonical, 'mean_event_rise19_sec', keyCols, 0);
cellCanonical.mean_event_fall91_sec_norm = localNormByCellBaseline(cellCanonical, 'mean_event_fall91_sec', keyCols, 0);
cellCanonical.mean_event_decay_tau_sec_norm = localNormByCellBaseline(cellCanonical, 'mean_event_decay_tau_sec', keyCols, 0);
cellCanonical.mean_event_signal_auc_dff_sec_norm = localNormByCellBaseline(cellCanonical, 'mean_event_signal_auc_dff_sec', keyCols, 0);
cellCanonical.mean_event_peak_amplitude_dff_norm = localNormByCellBaseline(cellCanonical, 'mean_event_peak_amplitude_dff', keyCols, 0);
cellCanonical.mean_event_max_propagation_speed_um_per_sec_norm = localNormByCellBaseline(cellCanonical, 'mean_event_max_propagation_speed_um_per_sec', keyCols, 0);
end

function [eventBoundary, cellBoundary] = computeBoundaryTables(canonicalOut, canonicalNorm, epochSpec, escapeThresholdFrac)
eventBoundary = table();
cellBoundary = table();

if ~(isfield(canonicalOut, 'subset') && isstruct(canonicalOut.subset) && ...
        isfield(canonicalOut.subset, 'events') && istable(canonicalOut.subset.events) && ...
        isfield(canonicalOut.subset, 'event_cell_map') && istable(canonicalOut.subset.event_cell_map))
    return
end
Eraw = canonicalOut.subset.events;
Mraw = canonicalOut.subset.event_cell_map;
needMap = {'event_uid','cell_uid','is_primary_cell','frac_outside_cell','outside_cell_area_npix','cell_overlap_frac_of_event'};
if ~all(ismember(needMap, Mraw.Properties.VariableNames))
    return
end
if ~ismember('epoch_idx', Eraw.Properties.VariableNames) || ~ismember('epoch_name', Eraw.Properties.VariableNames)
    A = assignEpochsToEvents(Eraw, epochSpec, 'method', 'auto');
    Eraw = A.events;
end
keepJoin = intersect({'event_uid','animal_id','slice_num','file_id','epoch_idx','epoch_name'}, Eraw.Properties.VariableNames);
Jraw = outerjoin(Mraw(:, needMap), Eraw(:, keepJoin), 'Keys', {'event_uid'}, 'Type', 'left', 'MergeKeys', true);
Jraw = Jraw(logical(Jraw.is_primary_cell) & Jraw.epoch_idx > 0, :);
if isempty(Jraw)
    return
end
Jraw.primary_cell_uid = string(Jraw.cell_uid);
Jraw.is_escape_gt_threshold = isfinite(Jraw.frac_outside_cell) & Jraw.frac_outside_cell >= double(escapeThresholdFrac);
eventBoundary = Jraw(:, {'event_uid','animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name', ...
    'frac_outside_cell','outside_cell_area_npix','cell_overlap_frac_of_event','is_escape_gt_threshold'});

if isfield(canonicalNorm, 'subset') && isstruct(canonicalNorm.subset) && ...
        isfield(canonicalNorm.subset, 'events') && istable(canonicalNorm.subset.events) && ...
        isfield(canonicalNorm.subset, 'event_cell_map') && istable(canonicalNorm.subset.event_cell_map)
    En = canonicalNorm.subset.events;
    Mn = canonicalNorm.subset.event_cell_map;
    if ~ismember('epoch_idx', En.Properties.VariableNames) || ~ismember('epoch_name', En.Properties.VariableNames)
        A = assignEpochsToEvents(En, epochSpec, 'method', 'auto');
        En = A.events;
    end
    if all(ismember(needMap, Mn.Properties.VariableNames))
        keepJoinN = intersect({'event_uid','epoch_idx','epoch_name'}, En.Properties.VariableNames);
        Jn = outerjoin(Mn(:, needMap), En(:, keepJoinN), 'Keys', {'event_uid'}, 'Type', 'left', 'MergeKeys', true);
        Jn = Jn(logical(Jn.is_primary_cell) & Jn.epoch_idx > 0, :);
        if ~isempty(Jn)
            Jn.primary_cell_uid = string(Jn.cell_uid);
            eventBoundary = addNormColumns(eventBoundary, Jn, {'frac_outside_cell'});
        end
    end
end
if ~ismember('frac_outside_cell_norm', eventBoundary.Properties.VariableNames)
    eventBoundary.frac_outside_cell_norm = localNormByCellBaseline(eventBoundary, 'frac_outside_cell', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
end
eventBoundary.escape_threshold_frac = repmat(double(escapeThresholdFrac), height(eventBoundary), 1);
eventBoundary = attachFileCalibration(eventBoundary, canonicalOut);
if ismember('spatial_res', eventBoundary.Properties.VariableNames)
    eventBoundary.outside_cell_area_um2 = double(eventBoundary.outside_cell_area_npix) .* (double(eventBoundary.spatial_res).^2);
    eventBoundary.outside_cell_area_um2_norm = localNormByCellBaseline(eventBoundary, 'outside_cell_area_um2', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
end

cellBoundary = summarizeCellMetricTable(eventBoundary, ...
    {'frac_outside_cell','mean_frac_outside_origin_cell'});
if isempty(cellBoundary)
    return
end
cellBoundary.frac_escape_events_gt_threshold = localGroupedFraction(eventBoundary, {'animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name'}, 'is_escape_gt_threshold');
cellBoundary.escape_threshold_frac = repmat(double(escapeThresholdFrac), height(cellBoundary), 1);
cellBoundary.mean_frac_outside_origin_cell_norm = localNormByCellBaseline(cellBoundary, 'mean_frac_outside_origin_cell', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellBoundary.frac_escape_events_gt_threshold_norm = localNormByCellBaseline(cellBoundary, 'frac_escape_events_gt_threshold', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
end

function [eventPartition, cellRecruitment, sliceRecruitment] = computePartitionTables(canonicalOut, epochSpec)
eventPartition = table();
cellRecruitment = table();
sliceRecruitment = table();

if ~(isfield(canonicalOut, 'subset') && isstruct(canonicalOut.subset) && ...
        isfield(canonicalOut.subset, 'events') && istable(canonicalOut.subset.events) && ...
        isfield(canonicalOut.subset, 'event_cell_map') && istable(canonicalOut.subset.event_cell_map))
    return
end

E = canonicalOut.subset.events;
M = canonicalOut.subset.event_cell_map;
needMap = {'event_uid','cell_uid','is_primary_cell','cell_first_overlap_frame','cell_last_overlap_frame', ...
    'cell_overlap_duration_frames','cell_peak_overlap_npix','cell_peak_overlap_frac_of_cell', ...
    'cell_peak_overlap_frac_of_event_frame','cell_entry_frame_cov_thresh', ...
    'cell_entry_delay_from_primary_frames','cell_entry_delay_from_primary_cov_thresh_frames'};
if ~all(ismember(needMap, M.Properties.VariableNames))
    return
end
if ~ismember('epoch_idx', E.Properties.VariableNames) || ~ismember('epoch_name', E.Properties.VariableNames)
    A = assignEpochsToEvents(E, epochSpec, 'method', 'auto');
    E = A.events;
end
joinCols = intersect({'event_uid','animal_id','slice_num','file_id','epoch_idx','epoch_name','primary_cell_uid'}, E.Properties.VariableNames);
if ~ismember('event_uid', joinCols)
    return
end
J = outerjoin(M(:, needMap), E(:, joinCols), 'Keys', {'event_uid'}, 'Type', 'left', 'MergeKeys', true);
if ~all(ismember({'file_id','epoch_idx','epoch_name'}, J.Properties.VariableNames))
    return
end
J = J(J.epoch_idx > 0, :);
if isempty(J)
    return
end
if ~ismember('primary_cell_uid', J.Properties.VariableNames)
    J.primary_cell_uid = repmat("", height(J), 1);
end
J.primary_cell_uid = string(J.primary_cell_uid);
J.partition_cell_uid = string(J.cell_uid);
J.is_primary_cell = logical(J.is_primary_cell);
for vn = {'cell_first_overlap_frame','cell_last_overlap_frame','cell_overlap_duration_frames', ...
        'cell_peak_overlap_npix','cell_peak_overlap_frac_of_cell','cell_peak_overlap_frac_of_event_frame', ...
        'cell_entry_frame_cov_thresh','cell_entry_delay_from_primary_frames','cell_entry_delay_from_primary_cov_thresh_frames'}
    if ismember(vn{1}, J.Properties.VariableNames)
        J.(vn{1}) = double(J.(vn{1}));
    end
end
if ismember('cell_entry_delay_from_primary_cov_thresh_frames', J.Properties.VariableNames)
    delayForRecruit = J.cell_entry_delay_from_primary_cov_thresh_frames;
else
    delayForRecruit = J.cell_entry_delay_from_primary_frames;
end
J.is_recruited_partition = ~J.is_primary_cell & isfinite(delayForRecruit) & delayForRecruit > 0;
eventPartition = J(:, intersect({'event_uid','animal_id','slice_num','file_id','primary_cell_uid','partition_cell_uid', ...
    'epoch_idx','epoch_name','is_primary_cell','is_recruited_partition','cell_first_overlap_frame','cell_last_overlap_frame', ...
    'cell_overlap_duration_frames','cell_peak_overlap_npix','cell_peak_overlap_frac_of_cell', ...
    'cell_peak_overlap_frac_of_event_frame','cell_entry_frame_cov_thresh', ...
    'cell_entry_delay_from_primary_frames','cell_entry_delay_from_primary_cov_thresh_frames'}, J.Properties.VariableNames));
eventPartition = attachFileCalibration(eventPartition, canonicalOut);
if ismember('frame_rate_s', eventPartition.Properties.VariableNames)
    eventPartition.cell_entry_delay_from_primary_sec = convertFramesToSecColumn(eventPartition, 'cell_entry_delay_from_primary_frames', 'frame_rate_s');
    eventPartition.cell_overlap_duration_sec = convertFramesToSecColumn(eventPartition, 'cell_overlap_duration_frames', 'frame_rate_s');
    eventPartition.cell_entry_delay_from_primary_sec_norm = localNormByCellBaseline(eventPartition, 'cell_entry_delay_from_primary_sec', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
    eventPartition.cell_overlap_duration_sec_norm = localNormByCellBaseline(eventPartition, 'cell_overlap_duration_sec', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
    if ismember('cell_entry_delay_from_primary_cov_thresh_frames', eventPartition.Properties.VariableNames)
        eventPartition.cell_entry_delay_from_primary_cov_thresh_sec = convertFramesToSecColumn(eventPartition, 'cell_entry_delay_from_primary_cov_thresh_frames', 'frame_rate_s');
        eventPartition.cell_entry_delay_from_primary_cov_thresh_sec_norm = localNormByCellBaseline(eventPartition, 'cell_entry_delay_from_primary_cov_thresh_sec', {'animal_id','slice_num','file_id','primary_cell_uid'}, 0);
    end
end
if ismember('cell_peak_overlap_frac_of_cell', eventPartition.Properties.VariableNames)
    eventPartition.cell_peak_overlap_frac_of_cell_norm = localNormByCellBaseline(eventPartition, 'cell_peak_overlap_frac_of_cell', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
end
if ismember('cell_peak_overlap_frac_of_event_frame', eventPartition.Properties.VariableNames)
    eventPartition.cell_peak_overlap_frac_of_event_frame_norm = localNormByCellBaseline(eventPartition, 'cell_peak_overlap_frac_of_event_frame', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
end

primaryEvents = unique(eventPartition(:, {'event_uid','animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name'}), 'rows', 'stable');
if isempty(primaryEvents)
    return
end
nonPrimary = eventPartition(eventPartition.is_recruited_partition, :);
if isempty(nonPrimary)
    primaryEvents.n_neighbor_partitions = zeros(height(primaryEvents),1);
    primaryEvents.has_neighbor_recruitment = false(height(primaryEvents),1);
    primaryEvents.mean_neighbor_entry_delay_frames = NaN(height(primaryEvents),1);
    primaryEvents.mean_neighbor_peak_overlap_frac_of_cell = NaN(height(primaryEvents),1);
else
    [G, grp] = localGroupRows(nonPrimary, {'event_uid'});
    nG = max([G;0]);
    perEvent = grp;
    perEvent.n_neighbor_partitions = zeros(nG,1);
    perEvent.mean_neighbor_entry_delay_frames = NaN(nG,1);
    perEvent.mean_neighbor_peak_overlap_frac_of_cell = NaN(nG,1);
    for k = 1:nG
        Tk = nonPrimary(G == k,:);
        perEvent.n_neighbor_partitions(k) = height(Tk);
        x = double(Tk.cell_entry_delay_from_primary_frames);
        x = x(isfinite(x));
        if ~isempty(x)
            perEvent.mean_neighbor_entry_delay_frames(k) = mean(x, 'omitnan');
        end
        x = double(Tk.cell_peak_overlap_frac_of_cell);
        x = x(isfinite(x));
        if ~isempty(x)
            perEvent.mean_neighbor_peak_overlap_frac_of_cell(k) = mean(x, 'omitnan');
        end
    end
    perEvent.has_neighbor_recruitment = perEvent.n_neighbor_partitions > 0;
    primaryEvents = outerjoin(primaryEvents, perEvent(:, {'event_uid','n_neighbor_partitions','has_neighbor_recruitment', ...
        'mean_neighbor_entry_delay_frames','mean_neighbor_peak_overlap_frac_of_cell'}), ...
        'Keys', {'event_uid'}, 'Type', 'left', 'MergeKeys', true);
    for vn = {'n_neighbor_partitions','has_neighbor_recruitment','mean_neighbor_entry_delay_frames','mean_neighbor_peak_overlap_frac_of_cell'}
        if ~ismember(vn{1}, primaryEvents.Properties.VariableNames)
            if strcmp(vn{1}, 'has_neighbor_recruitment')
                primaryEvents.(vn{1}) = false(height(primaryEvents),1);
            else
                primaryEvents.(vn{1}) = NaN(height(primaryEvents),1);
            end
        end
    end
end
primaryEvents.n_neighbor_partitions(~isfinite(primaryEvents.n_neighbor_partitions)) = 0;
primaryEvents.has_neighbor_recruitment = logical(primaryEvents.has_neighbor_recruitment);

[G, grp] = localGroupRows(primaryEvents, {'animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name'});
nG = max([G;0]);
if nG == 0
    return
end
cellRecruitment = grp;
cellRecruitment.n_primary_events = zeros(nG,1);
cellRecruitment.n_events_with_neighbor_recruitment = zeros(nG,1);
cellRecruitment.frac_events_with_neighbor_recruitment = NaN(nG,1);
cellRecruitment.mean_neighbor_entry_delay_frames = NaN(nG,1);
cellRecruitment.mean_neighbor_peak_overlap_frac_of_cell = NaN(nG,1);
cellRecruitment.mean_neighbor_partitions_per_event = NaN(nG,1);
for k = 1:nG
    Tk = primaryEvents(G == k,:);
    cellRecruitment.n_primary_events(k) = height(Tk);
    cellRecruitment.n_events_with_neighbor_recruitment(k) = nnz(Tk.has_neighbor_recruitment);
    if height(Tk) > 0
        cellRecruitment.frac_events_with_neighbor_recruitment(k) = mean(double(Tk.has_neighbor_recruitment), 'omitnan');
        x = double(Tk.mean_neighbor_entry_delay_frames);
        x = x(isfinite(x));
        if ~isempty(x)
            cellRecruitment.mean_neighbor_entry_delay_frames(k) = mean(x, 'omitnan');
        end
        x = double(Tk.mean_neighbor_peak_overlap_frac_of_cell);
        x = x(isfinite(x));
        if ~isempty(x)
            cellRecruitment.mean_neighbor_peak_overlap_frac_of_cell(k) = mean(x, 'omitnan');
        end
        x = double(Tk.n_neighbor_partitions);
        x = x(isfinite(x));
        if ~isempty(x)
            cellRecruitment.mean_neighbor_partitions_per_event(k) = mean(x, 'omitnan');
        end
    end
end
cellRecruitment.frac_events_with_neighbor_recruitment_norm = localNormByCellBaseline(cellRecruitment, ...
    'frac_events_with_neighbor_recruitment', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellRecruitment.mean_neighbor_entry_delay_frames_norm = localNormByCellBaseline(cellRecruitment, ...
    'mean_neighbor_entry_delay_frames', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellRecruitment.mean_neighbor_peak_overlap_frac_of_cell_norm = localNormByCellBaseline(cellRecruitment, ...
    'mean_neighbor_peak_overlap_frac_of_cell', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellRecruitment.mean_neighbor_partitions_per_event_norm = localNormByCellBaseline(cellRecruitment, ...
    'mean_neighbor_partitions_per_event', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);
cellRecruitment = attachFileCalibration(cellRecruitment, canonicalOut);
if ismember('frame_rate_s', cellRecruitment.Properties.VariableNames)
    cellRecruitment.mean_neighbor_entry_delay_sec = convertFramesToSecColumn(cellRecruitment, 'mean_neighbor_entry_delay_frames', 'frame_rate_s');
else
    cellRecruitment.mean_neighbor_entry_delay_sec = NaN(height(cellRecruitment),1);
end
cellRecruitment.mean_neighbor_entry_delay_sec_norm = localNormByCellBaseline(cellRecruitment, ...
    'mean_neighbor_entry_delay_sec', {'animal_id','slice_num','file_id','primary_cell_uid'}, 1e-6);

[sG, sGrp] = localGroupRows(primaryEvents, {'animal_id','slice_num','file_id','epoch_idx','epoch_name'});
nSG = max([sG;0]);
sliceRecruitment = sGrp;
sliceRecruitment.n_primary_events = zeros(nSG,1);
sliceRecruitment.n_primary_events_with_neighbor_recruitment = zeros(nSG,1);
sliceRecruitment.frac_primary_events_with_neighbor_recruitment = NaN(nSG,1);
sliceRecruitment.mean_neighbor_entry_delay_frames = NaN(nSG,1);
sliceRecruitment.mean_neighbor_peak_overlap_frac_of_cell = NaN(nSG,1);
sliceRecruitment.mean_neighbor_partitions_per_event = NaN(nSG,1);
for k = 1:nSG
    Tk = primaryEvents(sG == k,:);
    sliceRecruitment.n_primary_events(k) = height(Tk);
    sliceRecruitment.n_primary_events_with_neighbor_recruitment(k) = nnz(Tk.has_neighbor_recruitment);
    if height(Tk) > 0
        sliceRecruitment.frac_primary_events_with_neighbor_recruitment(k) = mean(double(Tk.has_neighbor_recruitment), 'omitnan');
        x = double(Tk.mean_neighbor_entry_delay_frames); x = x(isfinite(x));
        if ~isempty(x), sliceRecruitment.mean_neighbor_entry_delay_frames(k) = mean(x, 'omitnan'); end
        x = double(Tk.mean_neighbor_peak_overlap_frac_of_cell); x = x(isfinite(x));
        if ~isempty(x), sliceRecruitment.mean_neighbor_peak_overlap_frac_of_cell(k) = mean(x, 'omitnan'); end
        x = double(Tk.n_neighbor_partitions); x = x(isfinite(x));
        if ~isempty(x), sliceRecruitment.mean_neighbor_partitions_per_event(k) = mean(x, 'omitnan'); end
    end
end
sliceRecruitment.frac_primary_events_with_neighbor_recruitment_norm = localNormByCellBaseline(sliceRecruitment, ...
    'frac_primary_events_with_neighbor_recruitment', {'animal_id','slice_num','file_id'}, 1e-6);
sliceRecruitment.mean_neighbor_entry_delay_frames_norm = localNormByCellBaseline(sliceRecruitment, ...
    'mean_neighbor_entry_delay_frames', {'animal_id','slice_num','file_id'}, 1e-6);
sliceRecruitment.mean_neighbor_peak_overlap_frac_of_cell_norm = localNormByCellBaseline(sliceRecruitment, ...
    'mean_neighbor_peak_overlap_frac_of_cell', {'animal_id','slice_num','file_id'}, 1e-6);
sliceRecruitment.mean_neighbor_partitions_per_event_norm = localNormByCellBaseline(sliceRecruitment, ...
    'mean_neighbor_partitions_per_event', {'animal_id','slice_num','file_id'}, 1e-6);
sliceRecruitment = attachFileCalibration(sliceRecruitment, canonicalOut);
if ismember('frame_rate_s', sliceRecruitment.Properties.VariableNames)
    sliceRecruitment.mean_neighbor_entry_delay_sec = convertFramesToSecColumn(sliceRecruitment, 'mean_neighbor_entry_delay_frames', 'frame_rate_s');
else
    sliceRecruitment.mean_neighbor_entry_delay_sec = NaN(height(sliceRecruitment),1);
end
sliceRecruitment.mean_neighbor_entry_delay_sec_norm = localNormByCellBaseline(sliceRecruitment, ...
    'mean_neighbor_entry_delay_sec', {'animal_id','slice_num','file_id'}, 1e-6);
end

function T = legacySpreadLookup(ctx, thresholdPct)
T = table();
if ~isfield(ctx, 'BWA_AQ2') || ~isstruct(ctx.BWA_AQ2) || ~isfield(ctx.BWA_AQ2, 'fts1')
    return
end
fts = ctx.BWA_AQ2.fts1;
if ~isfield(fts, 'fileID')
    return
end
fileID = string(fts.fileID(:));
n = numel(fileID);
cov = ctx.associatedCoverageTraces;
if isempty(cov)
    return
end
if numel(cov) < n
    n = numel(cov);
    fileID = fileID(1:n);
end

if thresholdPct > 1
    thr = double(thresholdPct) / 100;
else
    thr = double(thresholdPct);
    thresholdPct = thr * 100;
end

localEventID = zeros(n,1);
u = unique(fileID, 'stable');
for i = 1:numel(u)
    idx = find(fileID == u(i));
    localEventID(idx) = (1:numel(idx))';
end

spreadCount = NaN(n,1);
for i = 1:n
    tr = cov{i};
    if isempty(tr) || ~isnumeric(tr)
        spreadCount(i) = NaN;
        continue
    end
    if isvector(tr)
        tr = tr(:)';
    end
    mx = max(tr, [], 1, 'omitnan');
    spreadCount(i) = sum(mx >= thr);
end

T = table(fileID, localEventID, spreadCount, repmat(double(thresholdPct), n, 1), ...
    'VariableNames', {'file_id','local_event_id','spread_count','spread_threshold_pct'});
end

function sliceNetwork = computeSliceNetworkTables(canonicalOut, canonicalNorm, epochSpec, coincidenceWindowSec, burstBinSec, burstMinEvents, corrBinSec)
sliceNetwork = table();

if ~(isfield(canonicalOut, 'subset') && isstruct(canonicalOut.subset) && ...
        isfield(canonicalOut.subset, 'events') && istable(canonicalOut.subset.events) && ...
        isfield(canonicalOut.subset, 'cells') && istable(canonicalOut.subset.cells) && ...
        isfield(canonicalOut.subset, 'files') && istable(canonicalOut.subset.files))
    return
end

E = canonicalOut.subset.events;
C = canonicalOut.subset.cells;
F = canonicalOut.subset.files;
if ~ismember('epoch_idx', E.Properties.VariableNames) || ~ismember('epoch_name', E.Properties.VariableNames)
    A = assignEpochsToEvents(E, epochSpec, 'method', 'auto');
    E = A.events;
end
if ~ismember('primary_cell_uid', E.Properties.VariableNames)
    E.primary_cell_uid = repmat("", height(E), 1);
end
E.primary_cell_uid = string(E.primary_cell_uid);
E = E(E.epoch_idx > 0 & E.primary_cell_uid ~= "", :);
if isempty(E)
    return
end

cellBase = C(:, intersect({'animal_id','slice_num','file_id','cell_uid'}, C.Properties.VariableNames));
[cellG, cellGrp] = localGroupRows(cellBase, {'animal_id','slice_num','file_id'});
nCellG = max([cellG;0]);
cellCounts = cellGrp;
cellCounts.n_total_cells = zeros(nCellG,1);
for kk = 1:nCellG
    cellCounts.n_total_cells(kk) = nnz(cellG == kk);
end

fileCols = intersect({'file_id','frame_rate_s','n_frames'}, F.Properties.VariableNames);
T = outerjoin(E, F(:, fileCols), 'Keys', {'file_id'}, 'Type', 'left', 'MergeKeys', true);
T = outerjoin(T, cellCounts, 'Keys', intersect({'animal_id','slice_num','file_id'}, T.Properties.VariableNames), 'Type', 'left', 'MergeKeys', true);
if ~ismember('n_total_cells', T.Properties.VariableNames)
    T.n_total_cells = NaN(height(T),1);
end
if ~ismember('frame_rate_s', T.Properties.VariableNames)
    T.frame_rate_s = NaN(height(T),1);
end
if ~ismember('n_frames', T.Properties.VariableNames)
    T.n_frames = NaN(height(T),1);
end

keyCols = {'animal_id','slice_num','file_id','epoch_idx','epoch_name'};
[G, grp] = localGroupRows(T, keyCols);
nG = max([G;0]);
if nG == 0
    return
end

sliceNetwork = grp;
sliceNetwork.n_total_cells = NaN(nG,1);
sliceNetwork.n_active_cells = zeros(nG,1);
sliceNetwork.n_events = zeros(nG,1);
sliceNetwork.frac_active_cells = NaN(nG,1);
sliceNetwork.coincidence_rate = NaN(nG,1);
sliceNetwork.burst_fraction = NaN(nG,1);
sliceNetwork.mean_pairwise_cell_correlation = NaN(nG,1);
sliceNetwork.coincidence_window_sec = repmat(double(coincidenceWindowSec), nG, 1);
sliceNetwork.burst_bin_sec = repmat(double(burstBinSec), nG, 1);
sliceNetwork.burst_min_events = repmat(double(burstMinEvents), nG, 1);
sliceNetwork.corr_bin_sec = repmat(double(corrBinSec), nG, 1);

for k = 1:nG
    Tk = T(G == k,:);
    if isempty(Tk)
        continue
    end
    nTotal = double(Tk.n_total_cells(1));
    if ~(isfinite(nTotal) && nTotal > 0)
        nTotal = numel(unique(string(Tk.primary_cell_uid)));
    end
    sliceNetwork.n_total_cells(k) = nTotal;
    sliceNetwork.n_events(k) = height(Tk);
    sliceNetwork.n_active_cells(k) = numel(unique(string(Tk.primary_cell_uid)));
    if isfinite(nTotal) && nTotal > 0
        sliceNetwork.frac_active_cells(k) = sliceNetwork.n_active_cells(k) / nTotal;
    end

    tSec = localEventStartSec(Tk);
    cids = string(Tk.primary_cell_uid);
    sliceNetwork.coincidence_rate(k) = localCoincidenceRate(tSec, cids, double(coincidenceWindowSec));

    frameRateS = double(Tk.frame_rate_s(1));
    nFrames = double(Tk.n_frames(1));
    [epochStartSec, epochEndSec] = localEpochWindowSec(epochSpec, double(Tk.epoch_idx(1)), frameRateS, nFrames, tSec);
    sliceNetwork.burst_fraction(k) = localBurstFraction(tSec, epochStartSec, epochEndSec, double(burstBinSec), double(burstMinEvents));
    sliceNetwork.mean_pairwise_cell_correlation(k) = localMeanPairwiseCellCorrelation(tSec, cids, epochStartSec, epochEndSec, double(corrBinSec));
end

sliceNetwork.frac_active_cells_norm = localNormByCellBaseline(sliceNetwork, ...
    'frac_active_cells', {'animal_id','slice_num','file_id'}, 1e-6);
sliceNetwork.coincidence_rate_norm = localNormByCellBaseline(sliceNetwork, ...
    'coincidence_rate', {'animal_id','slice_num','file_id'}, 1e-6);
sliceNetwork.burst_fraction_norm = localNormByCellBaseline(sliceNetwork, ...
    'burst_fraction', {'animal_id','slice_num','file_id'}, 1e-6);
sliceNetwork.mean_pairwise_cell_correlation_norm = localNormByCellBaseline(sliceNetwork, ...
    'mean_pairwise_cell_correlation', {'animal_id','slice_num','file_id'}, 1e-6);

if isfield(canonicalNorm, 'subset') && isstruct(canonicalNorm.subset) && ...
        isfield(canonicalNorm.subset, 'events') && istable(canonicalNorm.subset.events)
    % Summary metrics are normalized from baseline locally; no direct event-level normalized join required.
end
end

function tSec = localEventStartSec(T)
tSec = NaN(height(T),1);
if ismember('t0_sec', T.Properties.VariableNames)
    tSec = double(T.t0_sec);
end
needFrame = ~isfinite(tSec);
if any(needFrame) && all(ismember({'t0_frame','frame_rate_s'}, T.Properties.VariableNames))
    fr = double(T.frame_rate_s);
    t0f = double(T.t0_frame);
    ok = needFrame & isfinite(fr) & fr > 0 & isfinite(t0f);
    tSec(ok) = t0f(ok) .* fr(ok);
end
end

function rate = localCoincidenceRate(tSec, cellIDs, winSec)
rate = NaN;
if nargin < 3 || ~(isfinite(winSec) && winSec >= 0)
    return
end
ok = isfinite(tSec) & string(cellIDs) ~= "";
tSec = double(tSec(ok));
cellIDs = string(cellIDs(ok));
n = numel(tSec);
if n == 0
    return
elseif n == 1
    rate = 0;
    return
end
coincident = false(n,1);
for i = 1:n
    dt = abs(tSec - tSec(i));
    coincident(i) = any(dt <= winSec & dt > 0 & cellIDs ~= cellIDs(i));
end
rate = mean(double(coincident), 'omitnan');
end

function frac = localBurstFraction(tSec, startSec, endSec, binSec, minEvents)
frac = NaN;
if ~(isfinite(binSec) && binSec > 0 && isfinite(minEvents) && minEvents > 0)
    return
end
if ~(isfinite(startSec) && isfinite(endSec) && endSec > startSec)
    return
end
edges = startSec:binSec:endSec;
if numel(edges) < 2 || edges(end) < endSec
    edges = [edges endSec]; %#ok<AGROW>
end
if numel(edges) < 2
    return
end
tSec = double(tSec(isfinite(tSec)));
counts = histcounts(tSec, edges);
if isempty(counts)
    return
end
frac = mean(double(counts >= minEvents), 'omitnan');
end

function rho = localMeanPairwiseCellCorrelation(tSec, cellIDs, startSec, endSec, binSec)
rho = NaN;
if ~(isfinite(binSec) && binSec > 0)
    return
end
if ~(isfinite(startSec) && isfinite(endSec) && endSec > startSec)
    return
end
ok = isfinite(tSec) & string(cellIDs) ~= "";
tSec = double(tSec(ok));
cellIDs = string(cellIDs(ok));
if numel(tSec) < 2
    return
end
edges = startSec:binSec:endSec;
if numel(edges) < 2 || edges(end) < endSec
    edges = [edges endSec]; %#ok<AGROW>
end
if numel(edges) < 3
    return
end
cells = unique(cellIDs, 'stable');
if numel(cells) < 2
    return
end
X = zeros(numel(edges)-1, numel(cells));
for ii = 1:numel(cells)
    X(:,ii) = histcounts(tSec(cellIDs == cells(ii)), edges)';
end
R = corr(X, 'Rows', 'pairwise');
if isempty(R)
    return
end
mask = triu(true(size(R)), 1) & isfinite(R);
vals = R(mask);
vals = vals(isfinite(vals));
if ~isempty(vals)
    rho = mean(vals, 'omitnan');
end
end

function [startSec, endSec] = localEpochWindowSec(epochSpec, epochIdx, frameRateS, nFrames, tSecObserved)
startSec = NaN;
endSec = NaN;
if isfield(epochSpec, 'frameRange') && ~isempty(epochSpec.frameRange) && size(epochSpec.frameRange,1) >= epochIdx
    frng = double(epochSpec.frameRange(epochIdx,:));
    if numel(frng) == 2 && all(isfinite(frng)) && isfinite(frameRateS) && frameRateS > 0
        startSec = frng(1) * frameRateS;
        if isinf(frng(2)) && isfinite(nFrames)
            endSec = nFrames * frameRateS;
        else
            endSec = frng(2) * frameRateS;
        end
        return
    end
end
if isfield(epochSpec, 'timeRange') && ~isempty(epochSpec.timeRange) && size(epochSpec.timeRange,1) >= epochIdx
    trng = double(epochSpec.timeRange(epochIdx,:));
    if numel(trng) == 2 && all(isfinite(trng))
        startSec = trng(1);
        endSec = trng(2);
        return
    end
end
tSecObserved = double(tSecObserved(isfinite(tSecObserved)));
if ~isempty(tSecObserved)
    startSec = min(tSecObserved);
    endSec = max(tSecObserved);
    if endSec <= startSec
        endSec = startSec + max(frameRateS, 1);
    end
end
end

function [Eraw, Enorm] = canonicalEventTables(canonicalOut, canonicalNorm, epochSpec, requiredCols)
Eraw = table();
Enorm = table();
if ~(isfield(canonicalOut, 'subset') && isstruct(canonicalOut.subset) && isfield(canonicalOut.subset, 'events') && istable(canonicalOut.subset.events))
    return
end
Eraw = canonicalOut.subset.events;
if ~ismember('epoch_idx', Eraw.Properties.VariableNames) || ~ismember('epoch_name', Eraw.Properties.VariableNames)
    A = assignEpochsToEvents(Eraw, epochSpec, 'method', 'auto');
    Eraw = A.events;
end
if ~ismember('primary_cell_uid', Eraw.Properties.VariableNames)
    Eraw.primary_cell_uid = repmat("", height(Eraw), 1);
end
Eraw.primary_cell_uid = string(Eraw.primary_cell_uid);
Eraw = Eraw(Eraw.epoch_idx > 0 & string(Eraw.primary_cell_uid) ~= "", :);
if ~all(ismember(requiredCols, Eraw.Properties.VariableNames))
    Eraw = table();
    return
end

if isfield(canonicalNorm, 'subset') && isstruct(canonicalNorm.subset) && isfield(canonicalNorm.subset, 'events') && istable(canonicalNorm.subset.events)
    Enorm = canonicalNorm.subset.events;
    if ~ismember('epoch_idx', Enorm.Properties.VariableNames) || ~ismember('epoch_name', Enorm.Properties.VariableNames)
        A = assignEpochsToEvents(Enorm, epochSpec, 'method', 'auto');
        Enorm = A.events;
    end
    if ~ismember('primary_cell_uid', Enorm.Properties.VariableNames)
        Enorm.primary_cell_uid = repmat("", height(Enorm), 1);
    end
    Enorm.primary_cell_uid = string(Enorm.primary_cell_uid);
    Enorm = Enorm(Enorm.epoch_idx > 0 & string(Enorm.primary_cell_uid) ~= "", :);
end
end

function Tout = addNormColumns(Tout, Tnorm, vars)
keyCols = {'event_uid','primary_cell_uid','epoch_idx'};
keyCols = keyCols(ismember(keyCols, Tout.Properties.VariableNames) & ismember(keyCols, Tnorm.Properties.VariableNames));
if isempty(keyCols)
    return
end
for i = 1:numel(vars)
    v = vars{i};
    if ~ismember(v, Tnorm.Properties.VariableNames)
        continue
    end
    rhs = Tnorm(:, [keyCols {v}]);
    rhs.Properties.VariableNames{end} = [v '_norm'];
    Tout = outerjoin(Tout, rhs, 'Keys', keyCols, 'Type', 'left', 'MergeKeys', true);
end
end

function Tcell = summarizeCellMetricTable(Tevent, metricPairs)
Tcell = table();
if isempty(Tevent)
    return
end
keyCols = {'animal_id','slice_num','file_id','primary_cell_uid','epoch_idx','epoch_name'};
[G, grp] = localGroupRows(Tevent, keyCols);
nG = max([G;0]);
if nG == 0
    return
end
Tcell = grp;
for i = 1:size(metricPairs,1)
    inCol = metricPairs{i,1};
    outCol = metricPairs{i,2};
    vals = NaN(nG,1);
    for k = 1:nG
        x = double(Tevent.(inCol)(G == k));
        x = x(isfinite(x));
        if ~isempty(x)
            vals(k) = mean(x, 'omitnan');
        end
    end
    Tcell.(outCol) = vals;
end
end

function frac = localGroupedFraction(T, keyCols, boolCol)
[G, ~] = localGroupRows(T, keyCols);
nG = max([G;0]);
frac = NaN(nG,1);
for k = 1:nG
    x = logical(T.(boolCol)(G == k));
    if ~isempty(x)
        frac(k) = mean(double(x), 'omitnan');
    end
end
end

function T = attachFileCalibration(T, canonicalOut)
if isempty(T) || ~ismember('file_id', T.Properties.VariableNames)
    return
end
if ~(isfield(canonicalOut, 'subset') && isstruct(canonicalOut.subset) && ...
        isfield(canonicalOut.subset, 'files') && istable(canonicalOut.subset.files))
    return
end
F = canonicalOut.subset.files;
keep = intersect({'file_id','frame_rate_s','spatial_res'}, F.Properties.VariableNames);
if ~ismember('file_id', keep)
    return
end
rhs = unique(F(:, keep), 'rows', 'stable');
joinKeep = keep(~ismember(keep, T.Properties.VariableNames));
if isempty(joinKeep)
    return
end
rhs = rhs(:, ['file_id', joinKeep]);
T = outerjoin(T, rhs, 'Keys', {'file_id'}, 'Type', 'left', 'MergeKeys', true);
end

function vals = convertPxToUmColumn(T, pxCol, spatialResCol)
vals = NaN(height(T), 1);
if ~endsWith(string(pxCol), "_px")
    error('convertPxToUmColumn:ExpectedPxColumn', ...
        'Expected a raw pixel column ending in _px, got %s.', pxCol);
end
if ~(ismember(pxCol, T.Properties.VariableNames) && ismember(spatialResCol, T.Properties.VariableNames))
    return
end
x = double(T.(pxCol));
r = double(T.(spatialResCol));
ok = isfinite(x) & isfinite(r);
vals(ok) = x(ok) .* r(ok);
end

function vals = convertFramesToSecColumn(T, frameCol, frameRateCol)
vals = NaN(height(T), 1);
if ~endsWith(string(frameCol), "_frames")
    error('convertFramesToSecColumn:ExpectedFrameColumn', ...
        'Expected a raw frame column ending in _frames, got %s.', frameCol);
end
if ~(ismember(frameCol, T.Properties.VariableNames) && ismember(frameRateCol, T.Properties.VariableNames))
    return
end
x = double(T.(frameCol));
r = double(T.(frameRateCol));
ok = isfinite(x) & isfinite(r);
vals(ok) = x(ok) .* r(ok);
end

function normVals = localNormByCellBaseline(T, metricCol, keyCols, epsVal)
normVals = NaN(height(T), 1);
if ~ismember(metricCol, T.Properties.VariableNames) || ~ismember('epoch_idx', T.Properties.VariableNames)
    return
end
if nargin < 4 || isempty(epsVal)
    epsVal = 0;
end
keyCols = keyCols(ismember(keyCols, T.Properties.VariableNames));
if isempty(keyCols)
    return
end
transform = localMetricNormTransform();

[G, grp] = localGroupRows(T, keyCols); %#ok<ASGLU>
nG = max([G;0]);
for k = 1:nG
    idx = (G == k);
    Tk = T(idx, :);
    if isempty(Tk)
        continue
    end
    b = double(Tk.(metricCol)(Tk.epoch_idx == 1));
    b = b(isfinite(b));
    if isempty(b)
        continue
    end
    x = double(Tk.(metricCol));
    if strcmp(transform, 'log2_fold')
        x = x + epsVal;
        b0 = bwaComputeBaselineValue(x, Tk.epoch_idx == 1, transform);
        base = repmat(b0, size(x));
        y = bwaApplyTransform(x, base, transform);
    else
        b0 = bwaComputeBaselineValue(x, Tk.epoch_idx == 1, transform);
        base = repmat(b0, size(x));
        y = bwaApplyTransform(x, base, transform);
    end
    normVals(idx) = y;
end
end

function out = localMetricNormTransform(varargin)
persistent transform
if nargin > 0
    transform = lower(char(string(varargin{1})));
end
if isempty(transform)
    transform = 'log2_fold';
end
out = transform;
end

function T = emptyRegistry()
T = table(string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
    string.empty(0,1), string.empty(0,1), false(0,1), ...
    string.empty(0,1), string.empty(0,1), ...
    string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
    string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
    string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
    string.empty(0,1), string.empty(0,1), string.empty(0,1), ...
    string.empty(0,1), ...
    'VariableNames', {'metric_id','unit_level','source_table','value_column_raw', ...
    'value_column_norm','supports_norm','default_comparison','default_plot', ...
    'public_metric_id','family','native_unit','display_unit','aggregation_path','aggregation_rule','vocabulary_notes', ...
    'storage_status','storage_level','metric_class','default_slice_summary','alternate_slice_summary','storage_notes'});
end

function T = addRegistryRow(T, metricID, unitLevel, sourceTable, rawCol, normCol, supportsNorm, defaultCmp, defaultPlot)
row = table(string(metricID), string(unitLevel), string(sourceTable), ...
    string(rawCol), string(normCol), logical(supportsNorm), ...
    string(defaultCmp), string(defaultPlot), ...
    strings(1,1), strings(1,1), strings(1,1), strings(1,1), strings(1,1), strings(1,1), strings(1,1), ...
    strings(1,1), strings(1,1), strings(1,1), strings(1,1), strings(1,1), strings(1,1), ...
    'VariableNames', T.Properties.VariableNames);
T = [T; row];
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
elseif isnumeric(v) || islogical(v) || iscategorical(v)
    s = string(v(:));
else
    s = string(v);
    s = s(:);
end
end
