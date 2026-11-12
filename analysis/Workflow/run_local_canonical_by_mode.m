function [out, epochSpec, epochMeta, presetName] = run_local_canonical_by_mode(pathIn, mode, opts, qc, loaded, makePlots)
% Shared mode routing for local workflows.

if nargin < 6
    makePlots = false;
end

presetName = '';
mode = lower(char(string(mode)));

switch mode
    case 'preset'
        if isempty(string(getFieldOr(opts, 'presetName', '')))
            [plan, planMeta] = loadAnalysisPlan(pathIn);
        else
            [plan, planMeta] = loadAnalysisPlan(pathIn, 'presetName', char(string(opts.presetName)));
        end

        presetName = plan.presetName;
        epochPreset = '';
        if ~isempty(string(getFieldOr(opts, 'epochPreset', '')))
            epochPreset = char(string(opts.epochPreset));
        end
        if isfield(plan, 'epochPreset') && ~isempty(plan.epochPreset)
            if isempty(strtrim(epochPreset))
                epochPreset = char(string(plan.epochPreset));
            end
        end
        [epochSpec, epochMeta, filesOverride] = resolveEpochSpec(pathIn, epochPreset, loaded, opts, ...
            'presetName', epochPreset);

        subset = getPlanField(plan, 'subset', struct());
        metrics = getPlanField(plan, 'metrics', getFieldOr(opts, 'metrics', {}));
        metricGroupBy = getPlanField(plan, 'metricGroupBy', {'animal_id','slice_num','file_id'});
        freqGroupBy = getPlanField(plan, 'freqGroupBy', metricGroupBy);
        areaPercentiles = getPlanField(plan, 'areaPercentiles', [0 100]);
        dataLimitRules = resolveDataLimitRules(opts, plan);

        out = runCanonicalAnalysisPipelineWithQC(pathIn, epochSpec, ...
            'subset', subset, ...
            'metrics', metrics, ...
            'freqGroupBy', freqGroupBy, ...
            'metricGroupBy', metricGroupBy, ...
            'areaPercentiles', areaPercentiles, ...
            'dataLimitRules', dataLimitRules, ...
            'filesOverride', filesOverride, ...
            'dropZeroRates', logical(getFieldOr(opts, 'dropZeroRates', true)), ...
            'makePlots', logical(makePlots), ...
            'qc_files', qc.qc_files, 'qc_cells', qc.qc_cells, 'qc_events', qc.qc_events, ...
            'strictQC', logical(getFieldOr(opts, 'strictQC', false)));
        out.plan = plan;
        out.planMeta = planMeta;
        out.epochMeta = epochMeta;
        out.dataLimitRules = dataLimitRules;

    otherwise
        epochPreset = char(string(getFieldOr(opts, 'epochPreset', '')));
        [epochSpec, epochMeta, filesOverride] = resolveEpochSpec(pathIn, epochPreset, loaded, opts, ...
            'presetName', epochPreset, ...
            'sz_t', firstOrNaN(loaded.files, 'sz_t'), ...
            'frameRate_s', firstOrNaN(loaded.files, 'frame_rate_s'));

        subset = struct();
        metricGroupBy = {'animal_id','slice_num','file_id'};
        freqGroupBy = metricGroupBy;
        dataLimitRules = resolveDataLimitRules(opts, struct());

        switch mode
            case 'genotype'
                if ~isempty(getFieldOr(opts, 'groupValues', {}))
                    subset.genotypes = cellstr(string(opts.groupValues));
                end
                metricGroupBy = {'genotype','animal_id','slice_num','file_id'};
                freqGroupBy = metricGroupBy;
            case 'treatment'
                if ~isempty(getFieldOr(opts, 'groupValues', {}))
                    subset.treatments = cellstr(string(opts.groupValues));
                end
                metricGroupBy = {'treatment','animal_id','slice_num','file_id'};
                freqGroupBy = metricGroupBy;
            case 'epoch'
                % default grouping already set
        end

        out = runCanonicalAnalysisPipelineWithQC(pathIn, epochSpec, ...
            'subset', subset, ...
            'metrics', getFieldOr(opts, 'metrics', {'area_px2','duration_frames'}), ...
            'freqGroupBy', freqGroupBy, ...
            'metricGroupBy', metricGroupBy, ...
            'dataLimitRules', dataLimitRules, ...
            'filesOverride', filesOverride, ...
            'dropZeroRates', logical(getFieldOr(opts, 'dropZeroRates', true)), ...
            'makePlots', logical(makePlots), ...
            'qc_files', qc.qc_files, 'qc_cells', qc.qc_cells, 'qc_events', qc.qc_events, ...
            'strictQC', logical(getFieldOr(opts, 'strictQC', false)));
        out.mode = mode;
        out.epochMeta = epochMeta;
        out.subsetRequested = subset;
        out.dataLimitRules = dataLimitRules;
end
end

function [epochSpec, epochMeta, filesOverride] = resolveEpochSpec(pathIn, requestedPreset, loaded, opts, varargin)
manifestFields = {'recordingManifestPath','epochManifestPath','windowManifestPath'};
provided = false(size(manifestFields));
for ii = 1:numel(manifestFields)
    provided(ii) = strlength(string(getFieldOr(opts, manifestFields{ii}, ''))) > 0;
end
if any(provided) && ~all(provided)
    error('run_local_canonical_by_mode:IncompleteEpochManifests', ...
        'recordingManifestPath, epochManifestPath, and windowManifestPath must be supplied together.');
end

filesOverride = table();
if all(provided)
    [epochSpec, filesOverride, epochMeta] = buildRecordingRelativeEpochSpec(loaded.files, ...
        char(string(opts.recordingManifestPath)), char(string(opts.epochManifestPath)), ...
        char(string(opts.windowManifestPath)), ...
        'includeRoles', getFieldOr(opts, 'epochIncludeRoles', {'baseline','response'}), ...
        'analysisSpecPath', char(string(getFieldOr(opts, 'analysisSpecPath', ''))));
    epochMeta.source = 'recording_manifests';
    epochMeta.requestedPreset = char(string(requestedPreset));
else
    [epochSpec, epochMeta] = loadEpochSpecWithFallback(pathIn, requestedPreset, varargin{:});
end
end

function rules = resolveDataLimitRules(opts, plan)
rules = {};

if isstruct(plan) && isfield(plan, 'dataLimitRules') && iscell(plan.dataLimitRules) && ~isempty(plan.dataLimitRules)
    rules = plan.dataLimitRules;
end

if isfield(opts, 'dataLimitRules') && iscell(opts.dataLimitRules) && ~isempty(opts.dataLimitRules)
    rules = opts.dataLimitRules;
    return
end

dt = "";
pp = [];
if isstruct(plan)
    if isfield(plan, 'dataLimitType') && ~isempty(plan.dataLimitType), dt = string(plan.dataLimitType); end
    if isfield(plan, 'dataLimitPercentiles') && ~isempty(plan.dataLimitPercentiles), pp = plan.dataLimitPercentiles; end
end
if isfield(opts, 'dataLimitType') && ~isempty(opts.dataLimitType), dt = string(opts.dataLimitType); end
if isfield(opts, 'dataLimitPercentiles') && ~isempty(opts.dataLimitPercentiles), pp = opts.dataLimitPercentiles; end

if isempty(rules) && strlength(dt) > 0 && ~isempty(pp) && ~isNoOpPercentileLimit(pp)
    rules = AddLimitingRuleset({}, 'percentiles', char(dt), pp);
end
end

function tf = isNoOpPercentileLimit(pp)
pp = double(pp(:))';
tf = numel(pp) == 2 && pp(1) <= 0 && pp(2) >= 100;
end

function v = getFieldOr(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end

function v = getPlanField(plan, f, d)
if isstruct(plan) && isfield(plan, f) && ~isempty(plan.(f))
    v = plan.(f);
else
    v = d;
end
end

function v = firstOrNaN(T, varName)
v = NaN;
if istable(T) && ismember(varName, T.Properties.VariableNames) && ~isempty(T.(varName))
    v = T.(varName)(1);
end
end

function [epochSpec, epochMeta] = loadEpochSpecWithFallback(pathIn, requestedPreset, varargin)
requestedPreset = char(string(requestedPreset));
try
    [epochSpec, epochMeta] = loadEpochSpecFromConfig(pathIn, varargin{:});
catch ME
    isPresetMissing = contains(char(string(ME.identifier)), 'PresetNotFound') || ...
        contains(char(string(ME.message)), 'Preset "') || ...
        contains(char(string(ME.message)), 'preset');
    if ~isempty(strtrim(requestedPreset)) && isPresetMissing
        warning('run_local_canonical_by_mode:EpochPresetFallback', ...
            'Epoch preset "%s" not found in analysisConfig. Falling back to analysisConfig defaultPreset.', ...
            requestedPreset);
        args = varargin;
        args = stripPresetArg(args);
        if nargout > 1
            [epochSpec, epochMeta] = loadEpochSpecFromConfig(pathIn, args{:});
        else
            epochSpec = loadEpochSpecFromConfig(pathIn, args{:});
            epochMeta = struct();
        end
        epochMeta.requestedPreset = requestedPreset;
        epochMeta.fallbackUsed = true;
    else
        rethrow(ME);
    end
end
end

function out = stripPresetArg(args)
out = {};
i = 1;
while i <= numel(args)
    key = string(args{i});
    if isscalar(key) && strcmpi(char(key), 'presetName')
        i = i + 2;
    else
        out{end+1} = args{i}; %#ok<AGROW>
        i = i + 1;
    end
end
end
