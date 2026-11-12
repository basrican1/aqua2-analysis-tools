function result = load_local_analysis_result(source)
% Normalize either a struct or the path to an exported analysis_result.mat
if istable(source)
    error('load_local_analysis_result:BadSource', 'Table inputs not supported.');
end
if isstruct(source) && isfield(source, 'frequency')
    result = source;
    return
end
if isstruct(source) && isfield(source, 'canonical')
    result = source;
    return
end
if isstruct(source) && isfield(source, 'result')
    result = source.result;
    return
end
if ischar(source) || (isstring(source) && isscalar(source))
    path = char(source);
    if isfolder(path)
        candidate = fullfile(path, 'analysis_result.mat');
    else
        candidate = path;
    end
    if ~exist(candidate, 'file')
        error('load_local_analysis_result:MissingFile', 'Could not load analysis_result.mat from %s', candidate);
    end
    tmp = load(candidate, 'result');
    if ~isfield(tmp, 'result')
        error('load_local_analysis_result:BadFile', 'File %s does not contain variable ''result''.', candidate);
    end
    result = tmp.result;
    return
end
error('load_local_analysis_result:BadSource', 'Unsupported source type for analysis result.');
end
