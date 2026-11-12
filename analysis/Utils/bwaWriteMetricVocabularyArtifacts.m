function out = bwaWriteMetricVocabularyArtifacts(exportDir, result)
%BWAWRITEMETRICVOCABULARYARTIFACTS Write metric vocabulary sidecars.

metricVocabulary = table();
if isstruct(result) && isfield(result, 'metricVocabulary') && istable(result.metricVocabulary)
    metricVocabulary = result.metricVocabulary;
elseif isstruct(result) && isfield(result, 'metrics') && isstruct(result.metrics) && ...
        isfield(result.metrics, 'registry') && istable(result.metrics.registry)
    metricVocabulary = bwaEnrichMetricRegistry(result.metrics.registry);
end

matPath = fullfile(exportDir, 'metric_vocabulary.mat');
csvPath = fullfile(exportDir, 'metric_vocabulary.csv');
save(matPath, 'metricVocabulary', '-v7.3');
if ~isempty(metricVocabulary)
    writetable(metricVocabulary, csvPath);
end

out = struct();
out.table = metricVocabulary;
out.csvPath = csvPath;
out.matPath = matPath;
end
