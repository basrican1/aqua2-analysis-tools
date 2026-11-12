function info = listLegacyModules(varargin)
% List valid module names and key config fields for runLegacyComputeThenPlot.
%
% Usage:
%   info = listLegacyModules();
%   info = listLegacyModules('verbose', true);

p = inputParser;
p.StructExpand = false;
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});

mods = struct();

mods.standard = struct();
mods.standard.name = 'standard';
mods.standard.aliases = {'bwa_aqua2_standard'};
mods.standard.compute = {'prepareLegacyAqua2Context'};
mods.standard.plot = {'plotBWA_Aqua2StandardPanels'};
mods.standard.configFields = {'modules', 'doPlots', 'returnCompute', 'subset', 'useQC'};

mods.polars = struct();
mods.polars.name = 'polars';
mods.polars.aliases = {'bwa_polars', 'bwa_aqua2_polars'};
mods.polars.compute = {'computeBWA_Aqua2PolarsBySliceData'};
mods.polars.plot = {'plotBWA_Aqua2PolarsBySlice'};
mods.polars.configFields = {'includeStandardFeatures','spiderInputScale','spiderSymmetricLimit','exportCompanionTable','companionTableVarName'};

mods.baviolins = struct();
mods.baviolins.name = 'baviolins';
mods.baviolins.aliases = {'ba_violin_batch'};
mods.baviolins.compute = {'computeBAviolinData'};
mods.baviolins.plot = {'plotBAviolinPanels'};
mods.baviolins.configFields = {'baViolinDatTrees','baDoByCellOnly','baOptFigSuffix'};

globalCfg = { ...
    'runDir', ...
    'useQC', 'qc_files', 'qc_cells', 'qc_events', 'strictQC', ...
    'subset', 'reduceToSubsetFiles', ...
    'epochSpec', 'epochPreset', 'frameRange', 'timeRange', 'EpochList', ...
    'memberlist', 'associatedCoverageTraces', 'memberCoverage', ...
    'modules', 'doPlots', 'returnCompute', 'plotStyle', 'statsConfig'};

info = struct();
info.modules = mods;
info.globalConfigFields = globalCfg;

if p.Results.verbose
    fprintf('\nAvailable modules for runLegacyComputeThenPlot:\n');
    printModule(mods.standard);
    printModule(mods.polars);
    printModule(mods.baviolins);

    fprintf('\nGlobal cfg fields:\n');
    for i = 1:numel(globalCfg)
        fprintf('  - %s\n', globalCfg{i});
    end
    fprintf('\n');
end
end

function printModule(m)
fprintf('\n- %s\n', m.name);
fprintf('  aliases: %s\n', strjoin(m.aliases, ', '));
fprintf('  compute: %s\n', strjoin(m.compute, ', '));
fprintf('  plot: %s\n', strjoin(m.plot, ', '));
fprintf('  key cfg: %s\n', strjoin(m.configFields, ', '));
end
