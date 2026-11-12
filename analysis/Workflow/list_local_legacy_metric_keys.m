function out = list_local_legacy_metric_keys(source)
% List available legacy-panel metric names from a local analysis result.
%
% Usage:
%   out = list_local_legacy_metric_keys(result)
%   out = list_local_legacy_metric_keys('.../analysis_result.mat')
%
% Related family-based selection for legacy panels:
%   plot_local_analysis_legacy_panels(result, 'metricFamilies', {'shape'})
%   plot_local_analysis_legacy_panels(result, 'metricFamilies', {'propagation'})
%   plot_local_analysis_legacy_panels(result, 'metricFamilies', {'recruitment'})
%
% Common new cell metric keys:
%   cell.mean_perimeter_um
%   cell.mean_circularity
%   cell.frac_propagative_events
%   cell.mean_propagation_distance_um
%   cell.mean_frac_outside_origin_cell
%   cell.frac_escape_events_gt_threshold
%   cell.mean_soma_centroid_dist_um
%   cell.mean_soma_start_dist_um
%   cell.mean_soma_end_dist_um
%   cell.frac_events_with_neighbor_recruitment
%   cell.mean_neighbor_entry_delay_sec
%   cell.mean_neighbor_peak_overlap_frac_of_cell
%   cell.mean_neighbor_partitions_per_event
%
% Output struct fields:
%   out.metricListKeys   - keys accepted by plot_local_analysis_legacy_panels(..., 'metricList', ...)
%   out.extraKeys        - optional extra-derived panel keys
%   out.cellMetricKeys   - optional cell-level metric IDs
%   out.table            - combined listing table

R = load_local_analysis_result(source);

metricListKeys = {};
extraKeys = {};
cellMetricKeys = {};

if isfield(R, 'derived') && isstruct(R.derived) && ...
        isfield(R.derived, 'standard') && isstruct(R.derived.standard)
    S = R.derived.standard;
    if isfield(S, 'baviolin') && isstruct(S.baviolin)
        metricListKeys = fieldnames(S.baviolin)';
    end

    if isfield(S, 'footprint') && isstruct(S.footprint)
        extraKeys{end+1} = 'footprint'; %#ok<AGROW>
    end
    if isfield(S, 'propagation') && isstruct(S.propagation)
        if isfield(S.propagation, 'distance') && isstruct(S.propagation.distance)
            extraKeys{end+1} = 'propagation.distance'; %#ok<AGROW>
        end
        if isfield(S.propagation, 'angle') && isstruct(S.propagation.angle)
            extraKeys{end+1} = 'propagation.angle'; %#ok<AGROW>
        end
    end
    if isfield(S, 'centroid') && isstruct(S.centroid) && ...
            isfield(S.centroid, 'distComp') && isstruct(S.centroid.distComp)
        extraKeys{end+1} = 'centroid.distComp'; %#ok<AGROW>
    end
end

if isfield(R, 'metrics') && isstruct(R.metrics) && ...
        isfield(R.metrics, 'registry') && istable(R.metrics.registry)
    T = R.metrics.registry;
    if all(ismember({'metric_id','unit_level'}, T.Properties.VariableNames))
        m = strcmpi(string(T.unit_level), "cell");
        cellMetricKeys = cellstr(unique(string(T.metric_id(m)), 'stable'))';
        cellMetricKeys = filterPreferredCellMetricKeys(cellMetricKeys);
    end
end

% Build combined table.
allNames = [metricListKeys, extraKeys, cellMetricKeys];
allTypes = [repmat({'metricList'}, 1, numel(metricListKeys)), ...
            repmat({'extra'}, 1, numel(extraKeys)), ...
            repmat({'cell_metric'}, 1, numel(cellMetricKeys))];
if isempty(allNames)
    TT = table(strings(0,1), strings(0,1), 'VariableNames', {'name','type'});
else
    TT = table(string(allNames(:)), string(allTypes(:)), 'VariableNames', {'name','type'});
end

out = struct();
out.metricListKeys = metricListKeys;
out.extraKeys = extraKeys;
out.cellMetricKeys = cellMetricKeys;
out.table = TT;

% Print compact summary for interactive use.
fprintf('\nAvailable legacy-panel metric names:\n');
if isempty(metricListKeys)
    fprintf('  metricList keys: (none)\n');
else
    fprintf('  metricList keys (%d): %s\n', numel(metricListKeys), strjoin(metricListKeys, ', '));
end
if isempty(extraKeys)
    fprintf('  extra keys: (none)\n');
else
    fprintf('  extra keys (%d): %s\n', numel(extraKeys), strjoin(extraKeys, ', '));
end
if isempty(cellMetricKeys)
    fprintf('  cell metric keys: (none)\n');
else
    fprintf('  cell metric keys (%d): %s\n', numel(cellMetricKeys), strjoin(cellMetricKeys, ', '));
end
fprintf('\n');
end

function keys = filterPreferredCellMetricKeys(keys)
if isempty(keys)
    return
end
keys = string(keys(:));
pairs = [ ...
    "cell.mean_soma_centroid_dist_px", "cell.mean_soma_centroid_dist_um"; ...
    "cell.mean_soma_start_dist_px", "cell.mean_soma_start_dist_um"; ...
    "cell.mean_soma_end_dist_px", "cell.mean_soma_end_dist_um"; ...
    "cell.mean_perimeter_px", "cell.mean_perimeter_um"; ...
    "cell.mean_propagation_distance_px", "cell.mean_propagation_distance_um"; ...
    "cell.mean_neighbor_entry_delay_frames", "cell.mean_neighbor_entry_delay_sec"];
keep = true(size(keys));
for i = 1:size(pairs,1)
    raw = pairs(i,1); cal = pairs(i,2);
    if any(keys == cal)
        keep(keys == raw) = false;
    end
end
keys = cellstr(keys(keep))';
end
