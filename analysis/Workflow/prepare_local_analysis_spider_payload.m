function [result, polData, meta] = prepare_local_analysis_spider_payload(source, varargin)
% Prepare and optionally persist the spider payload used by spider plots.
%
% Usage:
%   result = prepare_local_analysis_spider_payload(result)
%   result = prepare_local_analysis_spider_payload(result, 'requireEntityValues', true)
%   [result, polData, meta] = prepare_local_analysis_spider_payload('.../analysis_result.mat', ...
%       'requireEntityValues', true, 'savePreparedResult', true)
%
% Name-Value:
%   'requireEntityValues' : true/false (default true). When true, ensure the
%                           payload includes pooled event/cell vectors needed
%                           for spider_mean comparisonLevel values like
%                           'by_event', 'by_cell', and 'by_branch'.
%   'savePreparedResult'  : true/false (default false). If source is a path to
%                           an exported analysis_result.mat (or its folder),
%                           overwrite that file with the prepared result so
%                           later spider calls do not rebuild the payload.

p = inputParser;
p.StructExpand = false;
addParameter(p, 'requireEntityValues', true, @islogical);
addParameter(p, 'savePreparedResult', false, @islogical);
parse(p, varargin{:});

[polData, result, meta] = ensure_local_polars_payload(source, ...
    'requireEntityValues', p.Results.requireEntityValues);

if p.Results.savePreparedResult
    resultPath = resolveAnalysisResultPath(source);
    if strlength(string(resultPath)) == 0
        error('prepare_local_analysis_spider_payload:NoWritablePath', ...
            ['savePreparedResult=true requires a path source (analysis_result.mat ' ...
             'or its export folder), not an in-memory struct.']);
    end
    save(char(resultPath), 'result', '-v7.3');
end
end

function resultPath = resolveAnalysisResultPath(source)
resultPath = "";
if ischar(source) || (isstring(source) && isscalar(source))
    candidate = char(source);
    if isfolder(candidate)
        candidate = fullfile(candidate, 'analysis_result.mat');
    end
    if exist(candidate, 'file')
        resultPath = string(candidate);
    end
end
end
