function stats = compute_local_mode_stats(out, mode, opts)
% Shared stats summary for local workflows.
% Frequency stats consume out.frequency as a compatibility view rebuilt from
% the authoritative result.metrics.tables.cell_frequency table.

stats = struct();
if ~isfield(out, 'frequency') || ~isfield(out.frequency, 'summaryByGroupEpoch')
    return
end

freqTbl = out.frequency.summaryByGroupEpoch;
compVar = resolveComparisonVar(mode, freqTbl);
if isempty(compVar)
    return
end

if ~isempty(string(getFieldOr(opts, 'epochForBetweenGroups', ''))) && ismember('epoch_name', freqTbl.Properties.VariableNames)
    keepEpoch = string(freqTbl.epoch_name) == string(opts.epochForBetweenGroups);
    freqTbl = freqTbl(keepEpoch, :);
end

statsCfg = struct();
if strcmp(mode, 'epoch') && logical(getFieldOr(opts, 'pairedEpochStats', true)) && ismember('file_id', freqTbl.Properties.VariableNames)
    statsCfg.paired = true;
    statsCfg.pairID = freqTbl.file_id;
end

if height(freqTbl) > 1
    stats.frequency = summarizeMetricStatsCanonical(freqTbl, 'mean_rate_events_per_min', compVar, 'statsConfig', statsCfg);
end

if ~isfield(out, 'metrics') || isempty(fieldnames(out.metrics))
    return
end
metricNames = fieldnames(out.metrics);
for i = 1:numel(metricNames)
    m = metricNames{i};
    Tm = out.metrics.(m).summaryByGroupEpoch;
    if ~isempty(string(getFieldOr(opts, 'epochForBetweenGroups', ''))) && ismember('epoch_name', Tm.Properties.VariableNames)
        keepEpoch = string(Tm.epoch_name) == string(opts.epochForBetweenGroups);
        Tm = Tm(keepEpoch, :);
    end

    valueVar = ['mean_' m];
    if ~ismember(valueVar, Tm.Properties.VariableNames)
        continue
    end
    cVar = resolveComparisonVar(mode, Tm);
    if isempty(cVar) || height(Tm) < 2
        continue
    end

    statsCfgM = struct();
    if strcmp(mode, 'epoch') && logical(getFieldOr(opts, 'pairedEpochStats', true)) && ismember('file_id', Tm.Properties.VariableNames)
        statsCfgM.paired = true;
        statsCfgM.pairID = Tm.file_id;
    end
    stats.metrics.(m) = summarizeMetricStatsCanonical(Tm, valueVar, cVar, 'statsConfig', statsCfgM);
end
end

function compVar = resolveComparisonVar(mode, T)
compVar = '';
switch lower(char(string(mode)))
    case 'epoch'
        if ismember('epoch_name', T.Properties.VariableNames), compVar = 'epoch_name'; end
    case 'genotype'
        if ismember('genotype', T.Properties.VariableNames), compVar = 'genotype'; end
    case 'treatment'
        if ismember('treatment', T.Properties.VariableNames), compVar = 'treatment';
        elseif ismember('drug', T.Properties.VariableNames), compVar = 'drug';
        elseif ismember('virus', T.Properties.VariableNames), compVar = 'virus';
        end
    case 'preset'
        if ismember('epoch_name', T.Properties.VariableNames), compVar = 'epoch_name';
        elseif ismember('genotype', T.Properties.VariableNames), compVar = 'genotype';
        elseif ismember('treatment', T.Properties.VariableNames), compVar = 'treatment';
        end
end
end

function v = getFieldOr(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end
