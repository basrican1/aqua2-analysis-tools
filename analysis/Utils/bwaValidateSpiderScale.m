function scale = bwaValidateSpiderScale(scaleIn, polData, msgPrefix)
%BWAVALIDATESPIDERSCALE Validate fold/log2 spider input scale.

scale = lower(char(string(scaleIn)));
if ~ismember(scale, {'fold','log2'})
    error('%s:BadScale', msgPrefix, 'spiderInputScale must be ''fold'' or ''log2''.');
end

if strcmp(scale, 'fold')
    v = [];
    if isfield(polData, 'valuesA') && isnumeric(polData.valuesA)
        v = [v; polData.valuesA(:)];
    end
    if isfield(polData, 'valuesC') && isnumeric(polData.valuesC)
        v = [v; polData.valuesC(:)];
    end
    v = v(isfinite(v));
    if ~isempty(v) && any(v <= 0)
        warning('%s:AutoScaleLog2', ...
            'Spider payload contains non-positive values; using spiderInputScale=''log2''.', msgPrefix);
        scale = 'log2';
    end
end
end
