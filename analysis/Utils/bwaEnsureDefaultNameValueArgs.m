function args = bwaEnsureDefaultNameValueArgs(argsIn, defaultPairs)
%BWAENSUREDEFAULTNAMEVALUEARGS Add missing default name-value pairs.

args = argsIn;
if ~iscell(args)
    args = {};
end
if ~iscell(defaultPairs) || isempty(defaultPairs)
    return
end

for i = 1:2:numel(defaultPairs)
    k = defaultPairs{i};
    v = defaultPairs{i+1};
    if ~bwaHasNameValueArg(args, k)
        args = [args {k, v}]; %#ok<AGROW>
    end
end
end
