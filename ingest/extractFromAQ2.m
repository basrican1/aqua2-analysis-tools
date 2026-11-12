function datOut = extractFromAQ2(fileDBInfo, IDStruct, fieldList)
% Fast legacy extraction for the fixed BWA_AQ2 schema used by interpretAQuA2.

    %#ok<INUSD>
    filepath = fileDBInfo.filepath;
    filename = fileDBInfo.filename;
    fileID = fileDBInfo.fileID;
    animalID = IDStruct.animalID;
    sliceNum = IDStruct.sliceNum;

    fullPath = strcat(filepath, filesep, filename);
    matFileObject = matfile(fullPath);
    resStruct = matFileObject.res;
    fts1 = resStruct.fts1;
    opts = resStruct.opts;

    incluLmk = getPathOr(fts1, 'region.cell.incluLmk', []);
    tempSize = size(incluLmk);
    if numel(tempSize) < 2
        tempSize(end+1:2) = 0;
    end
    nCells = tempSize(1);
    nSomas = tempSize(2);

    t0 = getPathOr(fts1, 'loc.t0', []);
    nEvents = numel(t0);

    datOut = struct();
    datOut.fts1.loc.t0 = t0;
    datOut.fts1.loc.t1 = getPathOr(fts1, 'loc.t1', []);
    datOut.fts1.loc.xSpa = getPathOr(fts1, 'loc.xSpa', []);
    datOut.fts1.loc.xSpaTemp = getPathOr(fts1, 'loc.xSpaTemp', []);

    datOut.fts1.basic.area = getPathOr(fts1, 'basic.area', []);

    datOut.fts1.curve.duration = getPathOr(fts1, 'curve.duration', []);
    datOut.fts1.curve.rise19 = getPathOr(fts1, 'curve.rise19', []);
    datOut.fts1.curve.fall91 = getPathOr(fts1, 'curve.fall91', []);
    datOut.fts1.curve.decayTau = getPathOr(fts1, 'curve.decayTau', []);
    datOut.fts1.curve.dffAUC = getPathOr(fts1, 'curve.dffAUC', []);
    datOut.fts1.curve.dffMax2 = getPathOr(fts1, 'curve.dffMax2', []);

    datOut.fts1.region.cell.memberIdx = getPathOr(fts1, 'region.cell.memberIdx', []);
    datOut.fts1.region.cell.mask = getPathOr(fts1, 'region.cell.mask', []);
    datOut.fts1.region.cell.incluLmk = incluLmk;
    datOut.fts1.region.landMark.mask = getPathOr(fts1, 'region.landMark.mask', {});
    datOut.fts1.region.landmarkDist.distAvg = getPathOr(fts1, 'region.landmarkDist.distAvg', []);
    datOut.fts1.region.landmarkDist.distMin = getPathOr(fts1, 'region.landmarkDist.distMin', []);
    datOut.fts1.region.landmarkDir.chgToward = getPathOr(fts1, 'region.landmarkDir.chgToward', []);
    datOut.fts1.region.landmarkDir.chgAway = getPathOr(fts1, 'region.landmarkDir.chgAway', []);
    datOut.fts1.region.landmarkDir.chgTowardBefReach = getPathOr(fts1, 'region.landmarkDir.chgTowardBefReach', []);
    datOut.fts1.region.landmarkDir.chgAwayAftReach = getPathOr(fts1, 'region.landmarkDir.chgAwayAftReach', []);

    datOut.fts1.propagation.maxPropSpeed = getPathOr(fts1, 'propagation.maxPropSpeed', []);
    datOut.fts1.propagation.propGrow = getPathOr(fts1, 'propagation.propGrow', {});
    datOut.fts1.propagation.propShrink = getPathOr(fts1, 'propagation.propShrink', {});

    datOut.opts.sz = getPathOr(opts, 'sz', []);
    datOut.opts.frameRate = getPathOr(opts, 'frameRate', []);
    datOut.opts.spatialRes = getPathOr(opts, 'spatialRes', []);

    datOut.fts1.fileID = repmat({fileID}, 1, nEvents);
    datOut.DB.fileID = {fileID};
    datOut.DB.animalID = {animalID};
    datOut.DB.sliceNum = {sliceNum};
    datOut.DB.nEvents = nEvents;
    datOut.DB.nCells = nCells;
    datOut.DB.nSomas = nSomas;
end

function out = getPathOr(s, dotted, defaultVal)
out = defaultVal;
cur = s;
parts = strsplit(dotted, '.');
for ii = 1:numel(parts)
    if ~isstruct(cur) || ~isfield(cur, parts{ii})
        return
    end
    cur = cur.(parts{ii});
end
out = cur;
end
