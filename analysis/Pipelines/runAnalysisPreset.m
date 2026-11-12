function out = runAnalysisPreset(pathIn, varargin)
% Compatibility wrapper: run one preset with full compute
% (canonical + derived) using run_local_analysis_compute.
%
% This keeps legacy call sites working while aligning behavior to the
% unified local workflow.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));

p = inputParser;
p.StructExpand = false;
addRequired(p, 'pathIn', @ischar);
addParameter(p, 'presetName', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'qc_files', table(), @(x) istable(x) || isempty(x));
addParameter(p, 'qc_cells', table(), @(x) istable(x) || isempty(x));
addParameter(p, 'qc_events', table(), @(x) istable(x) || isempty(x));
addParameter(p, 'useDefaultQCTables', true, @islogical);
parse(p, pathIn, varargin{:});

qc_files = p.Results.qc_files;
qc_cells = p.Results.qc_cells;
qc_events = p.Results.qc_events;
if ~p.Results.useDefaultQCTables && isempty(qc_files) && isempty(qc_cells) && isempty(qc_events)
    qc_files = table();
    qc_cells = table();
    qc_events = table();
end

result = run_local_analysis_compute( ...
    'pathIn', pathIn, ...
    'mode', 'preset', ...
    'presetName', char(string(p.Results.presetName)), ...
    'useQC', true, ...
    'qc_files', qc_files, ...
    'qc_cells', qc_cells, ...
    'qc_events', qc_events, ...
    'exportOutputs', false, ...
    'runStats', true);

% Keep old top-level canonical fields for compatibility.
out = result;
if isfield(result, 'canonical') && isstruct(result.canonical)
    c = result.canonical;
    if isfield(c, 'frequency'), out.frequency = c.frequency; end
    if isfield(c, 'metrics'), out.metrics = c.metrics; end
    if isfield(c, 'subset'), out.subset = c.subset; end
    if isfield(c, 'plan'), out.plan = c.plan; end
    if isfield(c, 'planMeta'), out.planMeta = c.planMeta; end
    if isfield(c, 'epochMeta'), out.epochMeta = c.epochMeta; end
end
end
