function row = bwaLookupMetricVocabulary(metricID, registry)
%BWALOOKUPMETRICVOCABULARY Return one-row vocabulary metadata table for a metric.

metricID = string(metricID);
row = table();

if nargin >= 2 && istable(registry) && ~isempty(registry) && ismember('metric_id', registry.Properties.VariableNames)
    reg = bwaEnrichMetricRegistry(registry);
    hit = find(string(reg.metric_id) == metricID, 1, 'first');
    if ~isempty(hit)
        cols = intersect({'metric_id','public_metric_id','family','native_unit','display_unit','aggregation_path','aggregation_rule','vocabulary_notes'}, ...
            reg.Properties.VariableNames, 'stable');
        row = reg(hit, cols);
        return
    end
end

cat = bwaMetricVocabularyCatalog(metricID);
if ~isempty(cat)
    row = cat(1,:);
end
end
