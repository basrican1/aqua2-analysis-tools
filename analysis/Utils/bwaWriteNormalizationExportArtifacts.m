function out = bwaWriteNormalizationExportArtifacts(exportDir, result)
%BWAWRITENORMALIZATIONEXPORTARTIFACTS Write normalization provenance sidecars.

[bundle, summaryTable, metricTable] = bwaBuildNormalizationExportBundle(result);

matPath = fullfile(exportDir, 'normalization_provenance.mat');
csvPath = fullfile(exportDir, 'normalization_provenance.csv');
metricCsvPath = fullfile(exportDir, 'normalization_metric_provenance.csv');

normalizationProvenance = bundle; %#ok<NASGU>
normalizationSummary = summaryTable; %#ok<NASGU>
normalizationMetrics = metricTable; %#ok<NASGU>
save(matPath, 'normalizationProvenance', 'normalizationSummary', 'normalizationMetrics', '-v7.3');
writetable(summaryTable, csvPath);
if ~isempty(metricTable)
    writetable(metricTable, metricCsvPath);
end

out = struct();
out.bundle = bundle;
out.summaryCsv = csvPath;
out.metricCsv = metricCsvPath;
out.matFile = matPath;
end
