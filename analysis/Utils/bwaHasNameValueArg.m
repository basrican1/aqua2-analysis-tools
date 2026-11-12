function tf = bwaHasNameValueArg(args, key)
%BWAHASNAMEVALUEARG True if key exists in name-value cell args.

tf = false;
if ~iscell(args) || isempty(args)
    return
end
for i = 1:2:numel(args)
    if i <= numel(args) && (ischar(args{i}) || isstring(args{i})) && strcmpi(string(args{i}), string(key))
        tf = true;
        return
    end
end
end
