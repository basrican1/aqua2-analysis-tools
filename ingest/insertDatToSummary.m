function outSummary = insertDatToSummary(inSummary, newDat, level)
    %#ok<INUSD>
    outSummary = inSummary;
    spec = legacyMergeSpec();
    for ii = 1:numel(spec)
        outSummary = mergeStructs(newDat, spec(ii).path, outSummary, spec(ii).path, spec(ii).dim);
    end
end

function spec = legacyMergeSpec()
paths = { ...
    'DB.fileID'
    'DB.animalID'
    'DB.sliceNum'
    'DB.nEvents'
    'DB.nCells'
    'DB.nSomas'
    'fts1.loc.t0'
    'fts1.loc.t1'
    'fts1.loc.xSpa'
    'fts1.loc.xSpaTemp'
    'fts1.basic.area'
    'fts1.curve.duration'
    'fts1.curve.rise19'
    'fts1.curve.fall91'
    'fts1.curve.decayTau'
    'fts1.curve.dffAUC'
    'fts1.curve.dffMax2'
    'fts1.region.cell.memberIdx'
    'fts1.region.cell.mask'
    'fts1.region.cell.incluLmk'
    'fts1.region.landMark.mask'
    'fts1.region.landmarkDist.distAvg'
    'fts1.region.landmarkDist.distMin'
    'fts1.region.landmarkDir.chgToward'
    'fts1.region.landmarkDir.chgAway'
    'fts1.region.landmarkDir.chgTowardBefReach'
    'fts1.region.landmarkDir.chgAwayAftReach'
    'fts1.propagation.maxPropSpeed'
    'fts1.propagation.propGrow'
    'fts1.propagation.propShrink'
    'fts1.fileID'
    'opts.sz'
    'opts.frameRate'
    'opts.spatialRes'};

spec = repmat(struct('path', '', 'dim', 2), numel(paths), 1);
for ii = 1:numel(paths)
    spec(ii).path = paths{ii};
    spec(ii).dim = legacyMergeDim(paths{ii});
end
end

function dim = legacyMergeDim(pathStr)
if contains(pathStr, 'memberIdx') || contains(pathStr, 'cell.mask') || ...
        contains(pathStr, 'landMark.mask') || contains(pathStr, 'incluLmk') || ...
        contains(pathStr, 'landmarkDist') || contains(pathStr, 'landmarkDir') || ...
        contains(pathStr, 'opts.sz')
    dim = 1;
else
    dim = 2;
end
end
