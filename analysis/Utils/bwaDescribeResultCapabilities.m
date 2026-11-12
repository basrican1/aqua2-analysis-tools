function caps = bwaDescribeResultCapabilities(result)
%BWADESCRIBERESULTCAPABILITIES Summarize which local-analysis payloads are available.

result = load_local_analysis_result(result);

caps = struct();
caps.hasCanonical = isfield(result, 'canonical') && isstruct(result.canonical);
caps.hasDerived = isfield(result, 'derived') && isstruct(result.derived);
caps.hasLegacyPayload = false;
caps.canRebuildFromPath = isfield(result, 'pathIn') && ~isempty(result.pathIn) && ...
    exist(char(string(result.pathIn)), 'dir') == 7;

if caps.hasDerived
    caps.hasDerivedContext = isfield(result.derived, 'context') && isstruct(result.derived.context);
    caps.hasDerivedStandard = isfield(result.derived, 'standard') && isstruct(result.derived.standard) && ...
        ~isempty(fieldnames(result.derived.standard));
    caps.hasPolarsPayload = isfield(result.derived, 'polars') && hasRequiredPolars(result.derived.polars);
else
    caps.hasDerivedContext = false;
    caps.hasDerivedStandard = false;
    caps.hasPolarsPayload = false;
end

if caps.hasDerivedContext && isfield(result.derived.context, 'BWA_AQ2') && ...
        isstruct(result.derived.context.BWA_AQ2) && ~isempty(fieldnames(result.derived.context.BWA_AQ2))
    caps.hasLegacyPayload = true;
elseif isfield(result, 'pathIn') && ~isempty(result.pathIn) && exist(char(string(result.pathIn)), 'dir') == 7
    try
        loaded = loadClusterOutputs(char(string(result.pathIn)), 'attachLegacyMeta', false, 'attachEventMeta', false);
        if isfield(loaded, 'hasLegacyPayload')
            caps.hasLegacyPayload = logical(loaded.hasLegacyPayload);
        elseif isfield(loaded, 'BWA_AQ2') && isstruct(loaded.BWA_AQ2) && ~isempty(fieldnames(loaded.BWA_AQ2))
            caps.hasLegacyPayload = true;
        end
    catch
    end
end

caps.canPlotSpider = caps.hasPolarsPayload || caps.hasLegacyPayload;
caps.canPlotLegacyPanels = caps.hasDerivedStandard;
caps.canPlotFileMetrics = caps.hasLegacyPayload;
caps.summary = buildSummary(caps);
end

function tf = hasRequiredPolars(polData)
need = {'labels','valuesA','valuesC','uniqueFiles'};
tf = isstruct(polData);
if ~tf
    return
end
for i = 1:numel(need)
    if ~isfield(polData, need{i}) || isempty(polData.(need{i}))
        tf = false;
        return
    end
end
end

function out = buildSummary(caps)
parts = strings(0,1);
if caps.hasCanonical
    parts(end+1,1) = "canonical"; %#ok<AGROW>
end
if caps.hasPolarsPayload
    parts(end+1,1) = "polars"; %#ok<AGROW>
end
if caps.hasDerivedStandard
    parts(end+1,1) = "derived_standard"; %#ok<AGROW>
end
if caps.hasLegacyPayload
    parts(end+1,1) = "legacy_payload"; %#ok<AGROW>
end
if isempty(parts)
    out = "none";
else
    out = strjoin(parts, ", ");
end
end
