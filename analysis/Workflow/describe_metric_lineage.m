function out = describe_metric_lineage(source, metricID)
%DESCRIBE_METRIC_LINEAGE Describe family and aggregation lineage for one or more metrics.

if nargin == 1
    metricID = source;
    source = [];
end

metricID = string(metricID(:));
out = table();

if ~isempty(source)
    result = load_local_analysis_result(source);
    if isstruct(result) && isfield(result, 'metrics') && isstruct(result.metrics) && ...
            isfield(result.metrics, 'registry') && istable(result.metrics.registry)
        reg = bwaEnrichMetricRegistry(result.metrics.registry);
        keep = ismember(string(reg.metric_id), metricID);
        if any(keep)
            cols = intersect({'metric_id','public_metric_id','family','native_unit','display_unit','aggregation_path','aggregation_rule', ...
                'storage_status','storage_level','metric_class','default_slice_summary','alternate_slice_summary','storage_notes','vocabulary_notes'}, ...
                reg.Properties.VariableNames, 'stable');
            out = reg(keep, cols);
        end
    end
end

if isempty(out) || ~ismember('metric_id', out.Properties.VariableNames)
    foundIDs = strings(0,1);
else
    foundIDs = string(out.metric_id);
end
missing = ~ismember(metricID, foundIDs);
if any(missing)
    cat = bwaMetricVocabularyCatalog(metricID(missing));
    if isempty(out)
        out = cat;
    elseif ~isempty(cat)
        add = table(string(cat.metric_id), string(cat.public_metric_id), string(cat.family), string(cat.native_unit), ...
            string(cat.display_unit), string(cat.aggregation_path), string(cat.aggregation_rule), ...
            string(cat.storage_status), string(cat.storage_level), string(cat.metric_class), ...
            string(cat.default_slice_summary), string(cat.alternate_slice_summary), string(cat.storage_notes), ...
            string(cat.notes), ...
            'VariableNames', out.Properties.VariableNames);
        out = [out; add];
    end
end

if nargout == 0
    disp(out);
end
end
