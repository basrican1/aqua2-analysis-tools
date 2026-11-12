function out = bwaWriteCapabilityExportArtifacts(exportDir, result)
%BWAWRITECAPABILITYEXPORTARTIFACTS Write result capability sidecars.

caps = bwaDescribeResultCapabilities(result);
summaryTable = capabilitySummaryTable(caps);

matPath = fullfile(exportDir, 'result_capabilities.mat');
csvPath = fullfile(exportDir, 'result_capabilities.csv');

resultCapabilities = caps; %#ok<NASGU>
resultCapabilitiesTable = summaryTable; %#ok<NASGU>
save(matPath, 'resultCapabilities', 'resultCapabilitiesTable', '-v7.3');
writetable(summaryTable, csvPath);

out = struct();
out.bundle = caps;
out.summaryCsv = csvPath;
out.matFile = matPath;
end

function T = capabilitySummaryTable(caps)
names = ["hasCanonical"; "hasLegacyPayload"; "hasPolarsPayload"; "hasDerivedStandard"; ...
    "canPlotSpider"; "canPlotLegacyPanels"; "canPlotFileMetrics"; "canRebuildFromPath"];
values = [logical(caps.hasCanonical); logical(caps.hasLegacyPayload); logical(caps.hasPolarsPayload); logical(caps.hasDerivedStandard); ...
    logical(caps.canPlotSpider); logical(caps.canPlotLegacyPanels); logical(caps.canPlotFileMetrics); logical(caps.canRebuildFromPath)];
notes = [ ...
    "Required base analysis contract"; ...
    "Legacy enrichment from BWA_AQ2.mat is available"; ...
    "Saved spider payload present"; ...
    "Legacy-derived standard payload present"; ...
    "Spider plotting available"; ...
    "Legacy BAviolin-style plotting available"; ...
    "File-level legacy plotting available"; ...
    "Result includes a usable pathIn for rebuild attempts" ...
    ];

T = table(names, values, notes, 'VariableNames', {'capability','available','notes'});
end
