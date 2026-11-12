function T = flatten_local_summary_significance(source, varargin)
% Flatten local summary significance output into one overview table.
%
% Examples:
%   S = compute_local_summary_significance(result, 'dataView', 'norm');
%   T = flatten_local_summary_significance(S);
%
%   T = flatten_local_summary_significance(result, ...
%       'dataView', 'norm', ...
%       'comparisonLevel', 'by_slice');

if isstruct(source) && isfield(source, 'source') && strcmp(source.source, 'compute_local_summary_significance')
    S = source;
else
    S = compute_local_summary_significance(source, varargin{:});
end

rows = table();
rows = [rows; collectRows('frequency', 'frequency', S.frequency)]; %#ok<AGROW>
rows = [rows; collectFamilyRows('canonical_metric', S.metrics)]; %#ok<AGROW>
rows = [rows; collectFamilyRows('derived_metric', S.derived)]; %#ok<AGROW>
rows = [rows; collectFamilyRows('cell_metric', S.cell_metrics, {'cell_frequency'})]; %#ok<AGROW>
rows = [rows; collectFamilyRows('slice_metric', S.slice_metrics)]; %#ok<AGROW>

if isempty(rows)
    T = table(strings(0,1), strings(0,1), strings(0,1), zeros(0,1), strings(0,1), ...
        zeros(0,1), zeros(0,1), false(0,1), strings(0,1), ...
        'VariableNames', {'family','metric','value_var','epoch_idx','epoch_name', ...
        'n','p_value','is_sig','method_used'});
    return
end

rows = filterPreferredSummaryRows(rows);
[~, ord] = sortrows(table(rows.family, rows.metric, rows.epoch_idx));
T = rows(ord, :);
end

function T = collectFamilyRows(familyName, S, excludeMetrics)
T = table();
if nargin < 3 || isempty(excludeMetrics)
    excludeMetrics = {};
end
if ~isstruct(S) || isempty(fieldnames(S))
    return
end
fn = fieldnames(S);
for i = 1:numel(fn)
    if any(strcmp(fn{i}, excludeMetrics))
        continue
    end
    T = [T; collectRows(familyName, fn{i}, S.(fn{i}))]; %#ok<AGROW>
end
end

function T = collectRows(familyName, metricName, item)
T = table();
if ~isstruct(item) || ~isfield(item, 'significance') || ~istable(item.significance)
    return
end
sig = item.significance;
if isempty(sig)
    return
end

n = zeros(height(sig),1);
p = NaN(height(sig),1);
isSig = false(height(sig),1);
method = repmat("", height(sig), 1);
valueVar = repmat("", height(sig), 1);
metric = repmat(string(metricName), height(sig), 1);
family = repmat(string(familyName), height(sig), 1);

if ismember('n', sig.Properties.VariableNames), n = double(sig.n); end
if ismember('p_value', sig.Properties.VariableNames), p = double(sig.p_value); end
if ismember('is_sig', sig.Properties.VariableNames), isSig = logical(sig.is_sig); end
if ismember('method_used', sig.Properties.VariableNames), method = string(sig.method_used); end
if isfield(item, 'valueVar')
    valueVar(:) = string(item.valueVar);
elseif isfield(item, 'metricCol')
    valueVar(:) = string(item.metricCol);
end

if ~ismember('epoch_name', sig.Properties.VariableNames)
    sig.epoch_name = "E" + string(sig.epoch_idx);
end

T = table(family, metric, valueVar, double(sig.epoch_idx), string(sig.epoch_name), ...
    n, p, isSig, method, 'VariableNames', ...
    {'family','metric','value_var','epoch_idx','epoch_name','n','p_value','is_sig','method_used'});
end

function T = filterPreferredSummaryRows(T)
if isempty(T)
    return
end
pairs = [ ...
    "area_px2", "area_um2"; ...
    "duration_frames", "duration_sec"; ...
    "rise19_frames", "rise19_sec"; ...
    "fall91_frames", "fall91_sec"; ...
    "decay_tau_frames", "decay_tau_sec"; ...
    "dff_auc", "dff_auc_sec"; ...
    "perimeter_px", "perimeter_um"; ...
    "propagation_distance_px", "propagation_distance_um"; ...
    "cell.mean_soma_centroid_dist_px", "cell.mean_soma_centroid_dist_um"; ...
    "cell.mean_soma_start_dist_px", "cell.mean_soma_start_dist_um"; ...
    "cell.mean_soma_end_dist_px", "cell.mean_soma_end_dist_um"; ...
    "cell.mean_perimeter_px", "cell.mean_perimeter_um"; ...
    "cell.mean_propagation_distance_px", "cell.mean_propagation_distance_um"; ...
    "cell.mean_neighbor_entry_delay_frames", "cell.mean_neighbor_entry_delay_sec"];
metric = string(T.metric);
keep = true(height(T),1);
for i = 1:size(pairs,1)
    raw = pairs(i,1); cal = pairs(i,2);
    if any(metric == cal)
        keep(metric == raw) = false;
    end
end
T = T(keep,:);
end
