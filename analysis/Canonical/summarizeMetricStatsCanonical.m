function out = summarizeMetricStatsCanonical(T, valueVar, groupVar, varargin)
% Run flexible group statistics on a canonical summary table.
%
% Example:
%   out = summarizeMetricStatsCanonical(metricTable, 'value', 'genotype', ...
%       'pairVar', 'animal_id', 'statsConfig', struct('paired', false));

p = inputParser;
p.StructExpand = false;
addRequired(p, 'T', @istable);
addRequired(p, 'valueVar', @(x) ischar(x) || isstring(x));
addRequired(p, 'groupVar', @(x) ischar(x) || isstring(x));
addParameter(p, 'pairVar', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'statsConfig', struct(), @isstruct);
parse(p, T, valueVar, groupVar, varargin{:});

valueVar = char(string(valueVar));
groupVar = char(string(groupVar));
pairVar = char(string(p.Results.pairVar));

if ~ismember(valueVar, T.Properties.VariableNames)
    error('summarizeMetricStatsCanonical:MissingValueVar', 'Missing value variable: %s', valueVar);
end
if ~ismember(groupVar, T.Properties.VariableNames)
    error('summarizeMetricStatsCanonical:MissingGroupVar', 'Missing group variable: %s', groupVar);
end

cfg = p.Results.statsConfig;
if ~isempty(pairVar)
    if ~ismember(pairVar, T.Properties.VariableNames)
        error('summarizeMetricStatsCanonical:MissingPairVar', 'Missing pair variable: %s', pairVar);
    end
    if ~isfield(cfg, 'pairID')
        cfg.pairID = T.(pairVar);
    end
end

values = T.(valueVar);
groups = T.(groupVar);

if isempty(fieldnames(cfg))
    stats = runGroupComparison(values, groups);
else
    nv = structToNV(cfg);
    stats = runGroupComparison(values, groups, nv{:});
end

out = struct();
out.valueVar = valueVar;
out.groupVar = groupVar;
out.pairVar = pairVar;
out.stats = stats;
out.summary = stats.summaryTable;
out.pairwise = stats.pairwiseC;
end

function nv = structToNV(s)
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
