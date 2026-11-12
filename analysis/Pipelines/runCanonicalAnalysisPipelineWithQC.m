function out = runCanonicalAnalysisPipelineWithQC(pathIn, epochSpec, varargin)
% Canonical-first analysis entrypoint with explicit QC application.
%
% This is equivalent to runCanonicalAnalysisPipeline, but inserts:
%   load canonical -> applyQCToCanonical -> subset -> compute -> optional plot

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

p = inputParser;
p.StructExpand = false;
addRequired(p, 'pathIn', @ischar);
addRequired(p, 'epochSpec', @isstruct);
addParameter(p, 'subset', struct(), @isstruct);
addParameter(p, 'metrics', {'area_px2','duration_frames','rise19_frames','fall91_frames','decay_tau_frames','dff_auc','dff_max2','max_prop_speed','perimeter_px','circularity','propagation_distance_px'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'freqGroupBy', {'file_id'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'metricGroupBy', {'animal_id','slice_num','file_id'}, @(x) iscell(x) || isstring(x));
addParameter(p, 'areaPercentiles', [0 100], @isnumeric);
addParameter(p, 'dataLimitRules', {}, @(x) iscell(x) || isempty(x));
addParameter(p, 'makePlots', false, @islogical);
addParameter(p, 'qc_files', table(), @(x) istable(x) || isempty(x));
addParameter(p, 'qc_cells', table(), @(x) istable(x) || isempty(x));
addParameter(p, 'qc_events', table(), @(x) istable(x) || isempty(x));
addParameter(p, 'strictQC', false, @islogical);
addParameter(p, 'filesOverride', table(), @(x) istable(x) || isempty(x));
addParameter(p, 'dropZeroRates', true, @islogical);
parse(p, pathIn, epochSpec, varargin{:});

loaded = loadClusterOutputs(pathIn, 'attachLegacyMeta', false, 'attachEventMeta', false);
if ~isempty(p.Results.filesOverride)
    loaded.files = applyFilesOverride(loaded.files, p.Results.filesOverride);
end

qcApplied = applyQCToCanonical(loaded.files, loaded.cells, loaded.events, loaded.event_cell_map, ...
    p.Results.qc_files, p.Results.qc_cells, p.Results.qc_events, 'strictQC', p.Results.strictQC);

subsetArgs = structToNameValue(p.Results.subset);
S = selectCanonicalSubset(qcApplied.files, qcApplied.cells, qcApplied.events, qcApplied.event_cell_map, subsetArgs{:});
S = applyCanonicalDataLimitRules(S, p.Results.dataLimitRules);

out = struct();
out.loaded = loaded;
out.qcApplied = qcApplied;
out.qcSummary = summarizeQCOutcome( ...
    struct('files', loaded.files, 'cells', loaded.cells, 'events', loaded.events, 'event_cell_map', loaded.event_cell_map), ...
    struct('files', qcApplied.files, 'cells', qcApplied.cells, 'events', qcApplied.events, 'event_cell_map', qcApplied.event_cell_map), ...
    p.Results.qc_files, p.Results.qc_cells, p.Results.qc_events);
out.subset = S;
out.dataLimitRules = p.Results.dataLimitRules;

freqAreaPct = p.Results.areaPercentiles;
if ~isempty(p.Results.dataLimitRules)
    % Global limits already applied to event table.
    freqAreaPct = [0 100];
end

out.frequency = computeNEventsPerMinCanonical(S.events, S.event_cell_map, S.files, epochSpec, ...
    'groupBy', p.Results.freqGroupBy, 'areaPercentiles', freqAreaPct, ...
    'dropZeroRates', p.Results.dropZeroRates, 'cells', S.cells);

metricNames = cellstr(string(p.Results.metrics));
out.metrics = struct();
for i = 1:numel(metricNames)
    m = metricNames{i};
    if ismember(m, S.events.Properties.VariableNames)
        out.metrics.(m) = computeMetricBySliceCanonical(S.events, S.files, m, epochSpec, ...
            'groupBy', p.Results.metricGroupBy);
    end
end

if p.Results.makePlots
    out.figures = struct();
    out.figures.frequency = plotNEventsPerMinCanonical(out.frequency);
    mNames = fieldnames(out.metrics);
    for i = 1:numel(mNames)
        m = mNames{i};
        out.figures.(m) = plotMetricBySliceCanonical(out.metrics.(m));
    end
end
end

function files = applyFilesOverride(canonicalFiles, overrideFiles)
if ~all(ismember({'file_id'}, canonicalFiles.Properties.VariableNames)) || ...
        ~all(ismember({'file_id'}, overrideFiles.Properties.VariableNames))
    error('runCanonicalAnalysisPipelineWithQC:FilesOverrideID', ...
        'Canonical and override file tables must contain file_id.');
end
if height(canonicalFiles) ~= height(overrideFiles) || ...
        ~all(ismember(string(canonicalFiles.file_id), string(overrideFiles.file_id)))
    error('runCanonicalAnalysisPipelineWithQC:FilesOverrideMismatch', ...
        'filesOverride must contain exactly the canonical file_id set.');
end
[matched, idx] = ismember(string(canonicalFiles.file_id), string(overrideFiles.file_id));
if ~all(matched)
    error('runCanonicalAnalysisPipelineWithQC:FilesOverrideMismatch', ...
        'filesOverride is missing one or more canonical file IDs.');
end
files = overrideFiles(idx,:);
end

function nv = structToNameValue(s)
if isempty(s)
    nv = {};
    return
end
f = fieldnames(s);
nv = cell(1, 2*numel(f));
for i = 1:numel(f)
    nv{2*i-1} = f{i};
    nv{2*i} = s.(f{i});
end
end

function S = applyCanonicalDataLimitRules(S, rules)
if isempty(rules) || ~iscell(rules) || isempty(S.events)
    return
end
E = S.events;
pass = true(height(E),1);
for i = 1:size(rules,1)
    if size(rules,2) < 3 || isempty(rules{i,3})
        continue
    end
    limitType = char(string(rules{i,1}));
    dataType = char(string(rules{i,2}));
    value = rules{i,3};
    col = mapDataTypeToEventCol(dataType);
    if isempty(col)
        error('runCanonicalAnalysisPipelineWithQC:UnknownDataLimitType', ...
            'Unknown canonical dataLimitType/dataType: %s. Use explicit canonical event columns such as area_px2 or duration_frames.', ...
            dataType);
    end
    if ~ismember(col, E.Properties.VariableNames)
        continue
    end
    x = double(E.(col));
    rulePass = applyRuleToVector(x, limitType, value);
    pass = pass & rulePass;
end

E = E(pass,:);
S.events = E;
if ismember('event_uid', E.Properties.VariableNames) && ismember('event_uid', S.event_cell_map.Properties.VariableNames)
    keepEvent = ismember(S.event_cell_map.event_uid, E.event_uid);
    S.event_cell_map = S.event_cell_map(keepEvent, :);
end
if ismember('cell_uid', S.event_cell_map.Properties.VariableNames) && ismember('cell_uid', S.cells.Properties.VariableNames)
    keepCell = ismember(S.cells.cell_uid, S.event_cell_map.cell_uid);
    S.cells = S.cells(keepCell, :);
end
end

function pass = applyRuleToVector(x, limitType, value)
pass = false(size(x));
finiteMask = isfinite(x);
if ~any(finiteMask)
    return
end
switch lower(limitType)
    case {'percentile','percentiles'}
        lo = prctile(x(finiteMask), value(1));
        hi = prctile(x(finiteMask), value(2));
    case {'threshold','thresholds'}
        lo = value(1);
        hi = value(2);
    otherwise
        lo = -inf;
        hi = inf;
end
pass = finiteMask & x >= lo & x <= hi;
end

function col = mapDataTypeToEventCol(dataType)
switch char(string(dataType))
    case 'area_px2'
        col = 'area_px2';
    case 'amplitude'
        col = 'dff_max2';
    case 'AUC'
        col = 'dff_auc';
    case 'duration_frames'
        col = 'duration_frames';
    case 'rise19_frames'
        col = 'rise19_frames';
    case 'fall91_frames'
        col = 'fall91_frames';
    case 'decay_tau_frames'
        col = 'decay_tau_frames';
    case 'maxspeed'
        col = 'max_prop_speed';
    otherwise
        col = '';
end
end
