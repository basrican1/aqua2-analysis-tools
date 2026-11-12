function exports = bwaWritePost2FrequencyExports(outputDir, result)
%BWAWRITEPOST2FREQUENCYEXPORTS Write tidy canonical frequency summaries.
% The cell table preserves zero-event cells. Recording and condition tables are
% descriptive summaries only; they do not perform inferential statistics.

exports = struct('cell_epoch_csv','', 'recording_epoch_csv','', 'condition_epoch_csv','');
if ~(isstruct(result) && isfield(result, 'canonical') && isstruct(result.canonical) && ...
        isfield(result.canonical, 'frequency') && isstruct(result.canonical.frequency) && ...
        isfield(result.canonical.frequency, 'ratesByCellEpoch') && ...
        istable(result.canonical.frequency.ratesByCellEpoch))
    return
end

T = result.canonical.frequency.ratesByCellEpoch;
if isempty(T) || ~ismember('rate_events_per_min', T.Properties.VariableNames)
    return
end
if exist(outputDir, 'dir') ~= 7, mkdir(outputDir); end

cellPath = fullfile(outputDir, 'post2_frequency_cell_epoch.csv');
writetable(T, cellPath);

recordingCols = presentColumns(T, {'file_id','recording_id','animal_id','slice_num', ...
    'bath_condition_id','stim_pattern','epoch_idx','epoch_name'});
conditionCols = presentColumns(T, {'bath_condition_id','stim_pattern','epoch_idx','epoch_name'});
recordingPath = fullfile(outputDir, 'post2_frequency_recording_epoch.csv');
conditionPath = fullfile(outputDir, 'post2_frequency_condition_epoch.csv');
writetable(summarizeRates(T, recordingCols), recordingPath);
writetable(summarizeRates(T, conditionCols), conditionPath);

exports.cell_epoch_csv = cellPath;
exports.recording_epoch_csv = recordingPath;
exports.condition_epoch_csv = conditionPath;
end

function cols = presentColumns(T, requested)
cols = requested(ismember(requested, T.Properties.VariableNames));
if isempty(cols)
    error('bwaWritePost2FrequencyExports:MissingGroups', 'No requested grouping columns are present.');
end
end

function S = summarizeRates(T, groupCols)
[G, keys] = findgroups(T(:, groupCols));
nGroups = height(keys);
S = keys;
S.n_cells = zeros(nGroups,1);
S.mean_rate_events_per_min = nan(nGroups,1);
S.median_rate_events_per_min = nan(nGroups,1);
S.std_rate_events_per_min = nan(nGroups,1);
S.sem_rate_events_per_min = nan(nGroups,1);
S.total_events = zeros(nGroups,1);
S.duration_min = nan(nGroups,1);
for ii = 1:nGroups
    rows = G == ii;
    x = double(T.rate_events_per_min(rows));
    x = x(isfinite(x));
    S.n_cells(ii) = numel(x);
    if ~isempty(x)
        S.mean_rate_events_per_min(ii) = mean(x);
        S.median_rate_events_per_min(ii) = median(x);
        S.std_rate_events_per_min(ii) = std(x);
        S.sem_rate_events_per_min(ii) = std(x) ./ sqrt(numel(x));
    end
    if ismember('n_events', T.Properties.VariableNames)
        S.total_events(ii) = sum(double(T.n_events(rows)), 'omitnan');
    end
    if ismember('duration_min', T.Properties.VariableNames)
        d = unique(double(T.duration_min(rows))); d = d(isfinite(d));
        if isscalar(d), S.duration_min(ii) = d; end
    end
end
end
