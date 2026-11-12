function BWA_sub = bwaReduceByFiles(BWA, fileIDs, filesTbl, eventsTbl)
%BWAREDUCEBYFILES Reduce legacy BWA struct to selected file IDs.
if isempty(fileIDs)
    BWA_sub = BWA;
    return
end
ids = cellstr(string(fileIDs(:)));
BWA_sub = ReduceAQ2ByFile(BWA, ids, 'filesTbl', filesTbl, 'eventsTbl', eventsTbl);
end
