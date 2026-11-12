function memberlist = bwaComputeMemberlistForBWA(BWA, opts, DB)
%BWACOMPUTEMEMBERLISTFORBWA Compute event-to-cell membership list.
[memberlist, ~, ~, ~] = eventBelongsCell(BWA.fts1, opts, DB);
end
