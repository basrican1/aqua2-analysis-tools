function outPath = writeLegacyPipelineConfigTemplate(pathIn, varargin)
% Write a starter config script for runLegacyComputeThenPlot().
%
% Example:
%   cfgPath = writeLegacyPipelineConfigTemplate(pathIn);
%   run(cfgPath);
%   out = runLegacyComputeThenPlot(pathIn, cfg);

p = inputParser;
p.StructExpand = false;
addRequired(p, 'pathIn', @ischar);
addParameter(p, 'filename', 'legacyPipelineConfig_template.m', @ischar);
parse(p, pathIn, varargin{:});

outPath = fullfile(pathIn, p.Results.filename);
fid = fopen(outPath, 'w');
if fid < 0
    error('writeLegacyPipelineConfigTemplate:OpenFailed', 'Could not open file for writing: %s', outPath);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '%% Template config for runLegacyComputeThenPlot\n');
fprintf(fid, 'cfg = struct();\n\n');
fprintf(fid, '%% Modules: standard | polars | baviolins\n');
fprintf(fid, 'cfg.modules = {''standard'', ''polars''};\n');
fprintf(fid, 'cfg.doPlots = true;\n');
fprintf(fid, 'cfg.returnCompute = true;\n\n');
fprintf(fid, '%% Optional global plot-style override (partial struct)\n');
fprintf(fid, 'cfg.plotStyle = struct();\n');
fprintf(fid, '%% Example:\n');
fprintf(fid, '%% s = getBWAPlotStyle();\n');
fprintf(fid, '%% s.fontSizeAxis = 11;\n');
fprintf(fid, '%% s.colors.groupA = [0.1 0.2 0.8];\n');
fprintf(fid, '%% cfg.plotStyle = s;\n');
fprintf(fid, '%% cfg.plotStyle = struct(''fontSizeAxis'', 11, ''lineWidthMain'', 2.0);\n\n');
fprintf(fid, '%% Optional statistics config used by BAviolin stats panels\n');
fprintf(fid, 'cfg.statsConfig = struct();\n');
fprintf(fid, '%% Examples:\n');
fprintf(fid, '%% cfg.statsConfig = struct(''method'',''auto'',''paired'',true,''parametric'',true);\n');
fprintf(fid, '%% cfg.statsConfig = struct(''method'',''ttest2'');\n\n');

fprintf(fid, '%% Optional QC\n');
fprintf(fid, 'cfg.useQC = false;\n');
fprintf(fid, 'cfg.qc_files = table();\n');
fprintf(fid, 'cfg.qc_cells = table();\n');
fprintf(fid, 'cfg.qc_events = table();\n');
fprintf(fid, 'cfg.strictQC = false;\n\n');

fprintf(fid, '%% Optional canonical subset (fields from selectCanonicalSubset)\n');
fprintf(fid, 'cfg.subset = struct();\n');
fprintf(fid, '%% Examples:\n');
fprintf(fid, '%% cfg.subset.genotypes = {''WT''};\n');
fprintf(fid, '%% cfg.subset.file_ids = {''file_001'', ''file_002''};\n');
fprintf(fid, 'cfg.reduceToSubsetFiles = true;\n\n');

fprintf(fid, '%% Epoch source: one of epochSpec or epochPreset or frame/time range\n');
fprintf(fid, 'cfg.epochSpec = struct();\n');
fprintf(fid, 'cfg.epochPreset = '''';\n');
fprintf(fid, 'cfg.frameRange = [];\n');
fprintf(fid, 'cfg.timeRange = [];\n');
fprintf(fid, 'cfg.EpochList = {};\n\n');

fprintf(fid, '%% Optional precomputed membership inputs (normally leave empty)\n');
fprintf(fid, 'cfg.memberlist = [];\n');
fprintf(fid, 'cfg.associatedCoverageTraces = {};\n');
fprintf(fid, 'cfg.memberCoverage = [];\n\n');

fprintf(fid, '%% Polars/spider options\n');
fprintf(fid, 'cfg.includeStandardFeatures = true;\n');
fprintf(fid, 'cfg.spiderInputScale = ''fold''; %% or ''log2''\n');
fprintf(fid, 'cfg.spiderSymmetricLimit = [];\n');
fprintf(fid, 'cfg.exportCompanionTable = true;\n');
fprintf(fid, 'cfg.companionTableVarName = ''BWA_spiderFoldTable'';\n\n');

fprintf(fid, '%% BAviolin batch options (used when modules includes ''baviolins'')\n');
fprintf(fid, 'cfg.baViolinDatTrees = {''basic.area'', ''curve.duration'', ''curve.rise19''};\n');
fprintf(fid, 'cfg.baDoByCellOnly = false;\n');
fprintf(fid, 'cfg.baOptFigSuffix = '''';\n\n');

fprintf(fid, '%% Optional explicit canonical run folder\n');
fprintf(fid, 'cfg.runDir = '''';\n');
end
