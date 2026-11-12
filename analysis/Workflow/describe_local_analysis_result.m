function out = describe_local_analysis_result(source)
%DESCRIBE_LOCAL_ANALYSIS_RESULT Summarize what a local-analysis result can support.

result = load_local_analysis_result(source);
caps = bwaDescribeResultCapabilities(result);

rows = {
    'canonical tables', logical(caps.hasCanonical), 'Required base analysis contract'
    'legacy payload', logical(caps.hasLegacyPayload), 'Needed for rebuilds from non-tabular AQuA2 data'
    'saved spider payload', logical(caps.hasPolarsPayload), 'Allows spider plotting without rebuild'
    'derived standard payload', logical(caps.hasDerivedStandard), 'Needed for legacy BAviolin-style panels'
    'can plot spider', logical(caps.canPlotSpider), 'Saved polars or legacy payload available'
    'can plot legacy panels', logical(caps.canPlotLegacyPanels), 'Requires derived standard payload'
    'can plot file metrics', logical(caps.canPlotFileMetrics), 'Requires legacy payload'
    };

T = cell2table(rows, 'VariableNames', {'capability','available','notes'});
out = struct();
out.capabilities = caps;
out.table = T;
out.summary = caps.summary;

if nargout == 0
    fprintf('Local analysis result capabilities\n');
    fprintf('  available payloads: %s\n', char(string(caps.summary)));
    for i = 1:height(T)
        fprintf('  %-22s %s\n', char(string(T.capability(i))), tfLabel(T.available(i)));
    end
end
end

function out = tfLabel(tf)
if tf
    out = 'yes';
else
    out = 'no';
end
end
