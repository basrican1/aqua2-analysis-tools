function out = runCanonicalAnalysisPipeline(pathIn, epochSpec, varargin)
% Canonical-first analysis entrypoint with compute/plot separation.
%
% Example:
%   loaded = loadClusterOutputs(pathIn, 'attachLegacyMeta', false, 'attachEventMeta', false);
%   ep = frameTimeConversionsByEpoch([1 300; 301 inf], [], loaded.files.sz_t(1), loaded.files.frame_rate_s(1), {'Base','Drug'});
%   out = runCanonicalAnalysisPipeline(pathIn, ep, ...
%       'subset', struct('genotypes', {{'WT'}}), ...
%       'metrics', {'area_px2','duration_frames','dff_max2'}, ...
%       'makePlots', true);

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
addParameter(p, 'makePlots', false, @islogical);
parse(p, pathIn, epochSpec, varargin{:});

loaded = loadClusterOutputs(pathIn, 'attachLegacyMeta', false, 'attachEventMeta', false);

subsetArgs = structToNameValue(p.Results.subset);
S = selectCanonicalSubset(loaded.files, loaded.cells, loaded.events, loaded.event_cell_map, subsetArgs{:});

out = struct();
out.loaded = loaded;
out.subset = S;

out.frequency = computeNEventsPerMinCanonical(S.events, S.event_cell_map, S.files, epochSpec, ...
    'groupBy', p.Results.freqGroupBy, 'areaPercentiles', p.Results.areaPercentiles);

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
