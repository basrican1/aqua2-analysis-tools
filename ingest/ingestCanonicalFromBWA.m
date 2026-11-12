function out = ingestCanonicalFromBWA(BWA_AQ2, DB_Files, DB_Animals, DB_ExpCond, varargin)
    % Build canonical tables from already-accumulated BWA_AQ2 (no raw file reread).

    BWA_setupPaths(fileparts(mfilename('fullpath')));

    p = inputParser;
    p.StructExpand = false;
    addRequired(p, 'BWA_AQ2', @isstruct);
    addRequired(p, 'DB_Files', @isstruct);
    addRequired(p, 'DB_Animals', @isstruct);
    addRequired(p, 'DB_ExpCond', @isstruct);
    addParameter(p, 'outputRoot', fullfile('derived', 'ingest_runs'), @ischar);
    addParameter(p, 'saveCsv', false, @islogical);
    addParameter(p, 'useParallel', false, @islogical);
    addParameter(p, 'minAssociatedCellCoverageFrac', 0.10, @isnumeric);
    parse(p, BWA_AQ2, DB_Files, DB_Animals, DB_ExpCond, varargin{:});

    outputRoot = p.Results.outputRoot;
    saveCsv = p.Results.saveCsv;
    useParallel = p.Results.useParallel;
    minAssociatedCellCoverageFrac = double(p.Results.minAssociatedCellCoverageFrac);

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    runDir = fullfile(outputRoot, timestamp);
    if ~exist(runDir, 'dir')
        mkdir(runDir);
    end

    filesOut = repmat(fileTemplate(), 0, 1);
    cellsOut = repmat(cellTemplate(), 0, 1);
    eventsOut = repmat(eventTemplate(), 0, 1);
    mapOut = repmat(mapTemplate(), 0, 1);

    if ~isfield(BWA_AQ2, 'DB') || isempty(BWA_AQ2.DB) || ~isfield(BWA_AQ2.DB, 'fileID')
        files = struct2table(filesOut);
        cells = struct2table(cellsOut);
        events = struct2table(eventsOut);
        event_cell_map = struct2table(mapOut);
    else
        bwaFileIDs = toCellstr(BWA_AQ2.DB.fileID);
        uniqueFileIDs = unique(bwaFileIDs, 'stable');

        allEventFileIDs = toCellstr(getFieldOr(BWA_AQ2.fts1, 'fileID', {}));
        all_t0 = forceCol(getFieldPathOr(BWA_AQ2.fts1, 'loc.t0', []));
        all_t1 = forceCol(getFieldPathOr(BWA_AQ2.fts1, 'loc.t1', []));
        all_area = forceCol(getFieldPathOr(BWA_AQ2.fts1, 'basic.area', []));
        all_duration = forceCol(getFieldPathOr(BWA_AQ2.fts1, 'curve.duration', []));
        all_rise19 = forceCol(getFieldPathOr(BWA_AQ2.fts1, 'curve.rise19', []));
        all_fall91 = forceCol(getFieldPathOr(BWA_AQ2.fts1, 'curve.fall91', []));
        all_decayTau = forceCol(getFieldPathOr(BWA_AQ2.fts1, 'curve.decayTau', []));
        all_dffAUC = forceCol(getFieldPathOr(BWA_AQ2.fts1, 'curve.dffAUC', []));
        all_dffMax2 = forceCol(getFieldPathOr(BWA_AQ2.fts1, 'curve.dffMax2', []));
        all_speed = forceCol(getFieldPathOr(BWA_AQ2.fts1, 'propagation.maxPropSpeed', []));
        all_propGrow = getFieldPathOr(BWA_AQ2.fts1, 'propagation.propGrow', {});
        all_propShrink = getFieldPathOr(BWA_AQ2.fts1, 'propagation.propShrink', {});
        all_memberIdx = getFieldPathOr(BWA_AQ2.fts1, 'region.cell.memberIdx', []);
        all_xSpa = getFieldPathOr(BWA_AQ2.fts1, 'loc.xSpa', {});
        all_xSpaTemp = getFieldPathOr(BWA_AQ2.fts1, 'loc.xSpaTemp', {});
        all_lmkDistAvg = getFieldPathOr(BWA_AQ2.fts1, 'region.landmarkDist.distAvg', []);
        all_lmkDistMin = getFieldPathOr(BWA_AQ2.fts1, 'region.landmarkDist.distMin', []);
        all_chgToward = getFieldPathOr(BWA_AQ2.fts1, 'region.landmarkDir.chgToward', []);
        all_chgAway = getFieldPathOr(BWA_AQ2.fts1, 'region.landmarkDir.chgAway', []);
        all_chgTowardBef = getFieldPathOr(BWA_AQ2.fts1, 'region.landmarkDir.chgTowardBefReach', []);
        all_chgAwayAft = getFieldPathOr(BWA_AQ2.fts1, 'region.landmarkDir.chgAwayAftReach', []);
        all_lmkMask = getFieldPathOr(BWA_AQ2.fts1, 'region.landMark.mask', {});
        all_incluLmk = getFieldPathOr(BWA_AQ2.fts1, 'region.cell.incluLmk', []);
        all_cellMask = getFieldPathOr(BWA_AQ2.fts1, 'region.cell.mask', {});
        all_cellCenter = getFieldPathOr(BWA_AQ2.fts1, 'region.cell.center', []);
        all_lmkCenter = getFieldPathOr(BWA_AQ2.fts1, 'region.landMark.center', []);
        nFiles = numel(uniqueFileIDs);
        cellOffsets = zeros(nFiles, 1);
        lmkOffsets = zeros(nFiles, 1);
        nCellsByFile = zeros(nFiles, 1);
        nSomasByFile = zeros(nFiles, 1);
        dbIdxByFile = zeros(nFiles, 1);
        cellOffset = 0;
        lmkOffset = 0;
        for ff = 1:nFiles
            fileID = uniqueFileIDs{ff};
            dbIdx = find(strcmp(bwaFileIDs, fileID), 1, 'first');
            dbIdxByFile(ff) = dbIdx;
            nCellsLocal = scalarCount(toScalarOrNaN(BWA_AQ2.DB.nCells, dbIdx));
            nSomasLocal = scalarCount(toScalarOrNaN(BWA_AQ2.DB.nSomas, dbIdx));
            cellOffsets(ff) = cellOffset;
            lmkOffsets(ff) = lmkOffset;
            nCellsByFile(ff) = nCellsLocal;
            nSomasByFile(ff) = nSomasLocal;
            cellOffset = cellOffset + nCellsLocal;
            lmkOffset = lmkOffset + nSomasLocal;
        end

        fileChunks = cell(nFiles, 1);
        cellChunks = cell(nFiles, 1);
        eventChunks = cell(nFiles, 1);
        mapChunks = cell(nFiles, 1);
        didParallel = false;
        if useParallel
            try
                if isempty(gcp('nocreate'))
                    parpool;
                end
                pObj = gcp('nocreate');
                disp(['Canonical: using parallel pool with ' num2str(pObj.NumWorkers) ' workers']);
                q = parallel.pool.DataQueue;
                nDone = 0;
                reportStep = max(1, floor(nFiles/20));
                afterEach(q, @canonicalProgressTick);
                parfor ff = 1:nFiles
                    [fileChunks{ff}, cellChunks{ff}, eventChunks{ff}, mapChunks{ff}] = ...
                        processOneBWAFile(ff, uniqueFileIDs, dbIdxByFile, cellOffsets, lmkOffsets, ...
                        nCellsByFile, nSomasByFile, BWA_AQ2, DB_Files, DB_Animals, DB_ExpCond, ...
                        allEventFileIDs, all_t0, all_t1, all_area, all_duration, all_rise19, all_fall91, ...
                        all_decayTau, all_dffAUC, all_dffMax2, all_speed, all_propGrow, all_propShrink, ...
                        all_memberIdx, all_xSpa, all_xSpaTemp, all_lmkDistAvg, all_lmkDistMin, ...
                        all_chgToward, all_chgAway, all_chgTowardBef, all_chgAwayAft, all_lmkMask, ...
                        all_incluLmk, all_cellMask, all_cellCenter, all_lmkCenter, minAssociatedCellCoverageFrac);
                    send(q, ff);
                end
                didParallel = true;
            catch ME
                disp('Canonical: parallel setup/execution failed');
                disp(ME.message);
                try
                    disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
                catch
                end
                didParallel = false;
            end
        end
        if ~didParallel
            if useParallel
                disp('Canonical: parallel unavailable, falling back to sequential mode');
            else
                disp('Canonical: using sequential mode');
            end
            for ff = 1:nFiles
                disp(['Canonical seq [' num2str(ff) '/' num2str(nFiles) ']']);
                [fileChunks{ff}, cellChunks{ff}, eventChunks{ff}, mapChunks{ff}] = ...
                    processOneBWAFile(ff, uniqueFileIDs, dbIdxByFile, cellOffsets, lmkOffsets, ...
                    nCellsByFile, nSomasByFile, BWA_AQ2, DB_Files, DB_Animals, DB_ExpCond, ...
                    allEventFileIDs, all_t0, all_t1, all_area, all_duration, all_rise19, all_fall91, ...
                    all_decayTau, all_dffAUC, all_dffMax2, all_speed, all_propGrow, all_propShrink, ...
                    all_memberIdx, all_xSpa, all_xSpaTemp, all_lmkDistAvg, all_lmkDistMin, ...
                    all_chgToward, all_chgAway, all_chgTowardBef, all_chgAwayAft, all_lmkMask, ...
                    all_incluLmk, all_cellMask, all_cellCenter, all_lmkCenter, minAssociatedCellCoverageFrac);
            end
        end
        for ff = 1:nFiles
            if ~isempty(fileChunks{ff}), filesOut = [filesOut; fileChunks{ff}]; end %#ok<AGROW>
            if ~isempty(cellChunks{ff}), cellsOut = [cellsOut; cellChunks{ff}]; end %#ok<AGROW>
            if ~isempty(eventChunks{ff}), eventsOut = [eventsOut; eventChunks{ff}]; end %#ok<AGROW>
            if ~isempty(mapChunks{ff}), mapOut = [mapOut; mapChunks{ff}]; end %#ok<AGROW>
        end

        files = struct2table(filesOut);
        cells = struct2table(cellsOut);
        events = struct2table(eventsOut);
        event_cell_map = struct2table(mapOut);
    end

    save(fullfile(runDir, 'files.mat'), 'files');
    save(fullfile(runDir, 'cells.mat'), 'cells');
    save(fullfile(runDir, 'events.mat'), 'events');
    save(fullfile(runDir, 'event_cell_map.mat'), 'event_cell_map');
    if saveCsv
        writetable(files, fullfile(runDir, 'files.csv'));
        writetable(cells, fullfile(runDir, 'cells.csv'));
        writetable(events, fullfile(runDir, 'events.csv'));
        writetable(event_cell_map, fullfile(runDir, 'event_cell_map.csv'));
    end

    manifest = struct();
    manifest.timestamp = timestamp;
    manifest.source = 'BWA_AQ2';
    manifest.run_dir = runDir;
    manifest.n_files = height(files);
    manifest.n_cells = height(cells);
    manifest.n_events = height(events);
    manifest.n_event_cell_map = height(event_cell_map);
    manifest.created_on = char(datetime('now'));
    save(fullfile(runDir, 'manifest.mat'), 'manifest');
    fid = fopen(fullfile(runDir, 'manifest.json'), 'w');
    if fid > 0
        fwrite(fid, jsonencode(manifest), 'char');
        fclose(fid);
    end
    fid = fopen(fullfile(outputRoot, 'latest_ingest.txt'), 'w');
    if fid > 0
        fprintf(fid, '%s\n', runDir);
        fclose(fid);
    end

    out = struct();
    out.runDir = runDir;
    out.files = files;
    out.cells = cells;
    out.events = events;
    out.event_cell_map = event_cell_map;
    out.manifest = manifest;

    function canonicalProgressTick(~)
        nDone = nDone + 1;
        if nDone == 1 || nDone == nFiles || mod(nDone, reportStep) == 0
            disp(['Canonical par [' num2str(nDone) '/' num2str(nFiles) ']']);
        end
    end
end

function [fileRows, cellRows, eventRows, mapRows] = processOneBWAFile(ff, uniqueFileIDs, dbIdxByFile, cellOffsets, lmkOffsets, ...
    nCellsByFile, nSomasByFile, BWA_AQ2, DB_Files, DB_Animals, DB_ExpCond, allEventFileIDs, all_t0, all_t1, all_area, ...
    all_duration, all_rise19, all_fall91, all_decayTau, all_dffAUC, all_dffMax2, all_speed, all_propGrow, all_propShrink, ...
    all_memberIdx, all_xSpa, all_xSpaTemp, all_lmkDistAvg, all_lmkDistMin, all_chgToward, all_chgAway, ...
    all_chgTowardBef, all_chgAwayAft, all_lmkMask, all_incluLmk, all_cellMask, all_cellCenter, all_lmkCenter, ...
    minAssociatedCellCoverageFrac)
fileRows = repmat(fileTemplate(), 0, 1);
cellRows = repmat(cellTemplate(), 0, 1);
eventRows = repmat(eventTemplate(), 0, 1);
mapRows = repmat(mapTemplate(), 0, 1);

fileID = uniqueFileIDs{ff};
dbIdx = dbIdxByFile(ff);
if isempty(dbIdx) || ~isfinite(dbIdx) || dbIdx < 1
    return
end
cellOffset = cellOffsets(ff);
lmkOffset = lmkOffsets(ff);
nCellsLocal = nCellsByFile(ff);
nSomasLocal = nSomasByFile(ff);

fileMeta = findFileMeta(DB_Files, fileID);
animalID = asChar(BWA_AQ2.DB.animalID{dbIdx});
animalMeta = findAnimalMeta(DB_Animals, animalID);
expMeta = findExpMeta(DB_ExpCond, fileID);
sliceNum = toScalarOrEmpty(BWA_AQ2.DB.sliceNum, dbIdx);
nEvents = toScalarOrNaN(BWA_AQ2.DB.nEvents, dbIdx);
nCells = toScalarOrNaN(BWA_AQ2.DB.nCells, dbIdx);
nSomas = toScalarOrNaN(BWA_AQ2.DB.nSomas, dbIdx);
localCellMask = sliceCellBlock(all_cellMask, cellOffset, nCellsLocal);
localLmkMask = sliceCellBlock(all_lmkMask, lmkOffset, nSomasLocal);
localIncluLmk = sliceNumericBlock(all_incluLmk, cellOffset, nCellsLocal, nSomasLocal);
localCellCenter = sliceNumericRows(all_cellCenter, cellOffset, nCellsLocal);
localLmkCenter = sliceNumericRows(all_lmkCenter, lmkOffset, nSomasLocal);
cellToSoma = NaN(0,1);

fr = toScalarOrNaN(BWA_AQ2.opts.frameRate, dbIdx);
spatialRes = toScalarOrNaN(BWA_AQ2.opts.spatialRes, dbIdx);
nFrames = NaN;
sz_x = NaN; sz_y = NaN; sz_z = NaN; sz_t = NaN;
if isfield(BWA_AQ2.opts, 'sz') && size(BWA_AQ2.opts.sz,1) >= dbIdx && size(BWA_AQ2.opts.sz,2) >= 4
    sz_x = BWA_AQ2.opts.sz(dbIdx,1);
    sz_y = BWA_AQ2.opts.sz(dbIdx,2);
    sz_z = BWA_AQ2.opts.sz(dbIdx,3);
    sz_t = BWA_AQ2.opts.sz(dbIdx,4);
    nFrames = BWA_AQ2.opts.sz(dbIdx,4);
end

frow = fileTemplate();
frow.file_id = string(fileID);
frow.animal_id = string(animalID);
frow.slice_num = string(toStr(sliceNum));
frow.filepath = string(fileMeta.filepath);
frow.filename = string(fileMeta.filename);
frow.exp_date = string(fileMeta.expDate);
frow.sz_x = sz_x;
frow.sz_y = sz_y;
frow.sz_z = sz_z;
frow.sz_t = sz_t;
frow.frame_rate_s = fr;
frow.n_frames = nFrames;
frow.spatial_res = spatialRes;
frow.n_events = nEvents;
frow.n_cells = nCells;
frow.n_somas = nSomas;
frow.animal_name = string(animalMeta.animalName);
frow.genotype = string(animalMeta.genotype);
frow.age = string(animalMeta.age);
frow.gender = string(animalMeta.gender);
frow.region_imaged = string(expMeta.regionImaged);
frow.virus = string(expMeta.virus);
frow.drug = string(expMeta.drug);
frow.stim = string(expMeta.stim);
frow.training = string(expMeta.training);
frow.image_epoc = string(expMeta.imageEpoc);
fileRows = frow;

if isnumeric(nCells) && ~isnan(nCells) && nCells > 0
    cellToSoma = resolveSingleSomaPerCellStandalone(localIncluLmk, localCellMask, localLmkMask, localCellCenter, localLmkCenter);
    cellRows = repmat(cellTemplate(), nCells, 1);
    for cc = 1:nCells
        crow = cellTemplate();
        crow.cell_uid = string(sprintf('%s_C%d', fileID, cc));
        crow.file_id = string(fileID);
        crow.animal_id = string(animalID);
        crow.slice_num = string(toStr(sliceNum));
        crow.local_cell_id = cc;
        crow.is_soma_like = false;
        crow.soma_landmark_idx = NaN;
        if cc <= numel(cellToSoma)
            crow.soma_landmark_idx = cellToSoma(cc);
            crow.is_soma_like = isfinite(crow.soma_landmark_idx) && crow.soma_landmark_idx > 0;
        end
        crow.cell_mask_npix = NaN;
        cellRows(cc,1) = crow;
    end
end

eventMask = strcmp(allEventFileIDs, fileID);
eventIdxAbs = find(eventMask);
nEventLocal = numel(eventIdxAbs);
evtPixCache = cell(nEventLocal, 1);
propGrowCache = cell(nEventLocal, 1);
propShrinkCache = cell(nEventLocal, 1);
xSpaTempCache = cell(nEventLocal, 1);
partitionFrameCache = cell(nEventLocal, 1);
for ee = 1:nEventLocal
    idx = eventIdxAbs(ee);
    evtPixCache{ee} = getEventSpatialPixels(all_xSpa, all_xSpaTemp, idx, BWA_AQ2.opts, dbIdx);
    propGrowCache{ee} = safeCellIndex(all_propGrow, idx);
    propShrinkCache{ee} = safeCellIndex(all_propShrink, idx);
    xSpaTempCache{ee} = safeCellIndex(all_xSpaTemp, idx);
    partitionFrameCache{ee} = canonicalEventFrameCache(xSpaTempCache{ee}, all_t0(idx), all_t1(idx), [sz_x sz_y sz_z sz_t]);
end
cellPixCache = cacheMaskPixels(localCellMask);
lmkPixCache = cacheMaskPixels(localLmkMask);
for ee = 1:nEventLocal
    idx = eventIdxAbs(ee);
    erow = eventTemplate();
    erow.event_uid = string(sprintf('%s_E%d', fileID, ee));
    erow.file_id = string(fileID);
    erow.animal_id = string(animalID);
    erow.slice_num = string(toStr(sliceNum));
    erow.local_event_id = ee;
    erow.t0_frame = safeIndex(all_t0, idx);
    erow.t1_frame = safeIndex(all_t1, idx);
    erow.t0_sec = erow.t0_frame * fr;
    erow.t1_sec = erow.t1_frame * fr;
    erow.area_px2 = safeIndex(all_area, idx);
    erow.duration_frames = safeIndex(all_duration, idx);
    erow.rise19_frames = safeIndex(all_rise19, idx);
    erow.fall91_frames = safeIndex(all_fall91, idx);
    erow.decay_tau_frames = safeIndex(all_decayTau, idx);
    erow.dff_auc = safeIndex(all_dffAUC, idx);
    erow.dff_max2 = safeIndex(all_dffMax2, idx);
    erow.max_prop_speed = safeIndex(all_speed, idx);
    evtPix = evtPixCache{ee};
    evMetric = canonicalEventMetrics(evtPix, [sz_x sz_y sz_z sz_t], [], propGrowCache{ee}, propShrinkCache{ee});
    erow.perimeter_px = evMetric.perimeter_px;
    erow.circularity = evMetric.circularity;
    erow.propagation_distance_px = evMetric.propagation_distance_px;
    erow.propagation_angle_deg = evMetric.propagation_angle_deg;

    [primaryCell, mappedCells, mappedScores] = extractMemberMapping(all_memberIdx, idx);
    owner = canonicalBranchOriginOwnership(BWA_AQ2, idx, 'candidateCells', mappedCells, 'cellMasks', localCellMask);
    branchOwnerCell = owner.ownerCell;
    branchXY = owner.originXY;
    branchFrame = owner.originFrame;
    if isfinite(branchOwnerCell) && branchOwnerCell > 0
        primaryCell = branchOwnerCell;
    else
        origin = canonicalOriginCell(xSpaTempCache{ee}, erow.t0_frame, erow.t1_frame, [sz_x sz_y sz_z sz_t], mappedCells, localCellMask);
        if isfinite(origin.originCell) && origin.originCell > 0
            primaryCell = origin.originCell;
        end
    end
    [mappedCells, mappedScores] = localMergePrimaryIntoCandidates(mappedCells, mappedScores, primaryCell);
    erow.branch_origin_x_px = branchXY(1);
    erow.branch_origin_y_px = branchXY(2);
    erow.branch_origin_frame = branchFrame;
    if isfinite(branchOwnerCell) && branchOwnerCell > 0
        erow.branch_origin_cell_uid = string(sprintf('%s_C%d', fileID, branchOwnerCell));
    else
        erow.branch_origin_cell_uid = "";
    end
    if ~isnan(primaryCell) && primaryCell > 0
        erow.primary_cell_uid = string(sprintf('%s_C%d', fileID, primaryCell));
    else
        erow.primary_cell_uid = "";
    end
    eventRows(end+1,1) = erow; %#ok<AGROW>

    eventMapRows = repmat(mapTemplate(), numel(mappedCells), 1);
    for mm = 1:numel(mappedCells)
        mrow = mapTemplate();
        mrow.event_uid = erow.event_uid;
        mrow.cell_uid = string(sprintf('%s_C%d', fileID, mappedCells(mm)));
        if numel(mappedScores) >= mm
            mrow.association_score = mappedScores(mm);
        else
            mrow.association_score = NaN;
        end
        mrow.is_primary_cell = (~isnan(primaryCell) && mappedCells(mm) == primaryCell);
        lmkIdx = NaN;
        if mappedCells(mm) <= numel(cellToSoma)
            lmkIdx = cellToSoma(mappedCells(mm));
        end
        mrow.soma_landmark_idx = lmkIdx;
        mrow.event_spatial_npix = numel(evtPix);
        cellPix = [];
        if numel(cellPixCache) >= mappedCells(mm)
            cellPix = cellPixCache{mappedCells(mm)};
        end
        mapMetric = canonicalEventMetrics(evtPix, [sz_x sz_y sz_z sz_t], cellPix, [], []);
        mrow.cell_overlap_npix = mapMetric.cell_overlap_npix;
        mrow.cell_overlap_frac_of_event = mapMetric.cell_overlap_frac_of_event;
        mrow.cell_overlap_frac_of_cell = mapMetric.cell_overlap_frac_of_cell;
        mrow.frac_outside_cell = mapMetric.frac_outside_cell;
        mrow.outside_cell_area_npix = mapMetric.outside_cell_area_npix;
        partMetric = canonicalEventCellPartitionMetrics(partitionFrameCache{ee}, all_t0(idx), all_t1(idx), [sz_x sz_y sz_z sz_t], cellPix, minAssociatedCellCoverageFrac);
        mrow.cell_first_overlap_frame = partMetric.cell_first_overlap_frame;
        mrow.cell_entry_frame_cov_thresh = partMetric.cell_entry_frame_cov_thresh;
        mrow.cell_last_overlap_frame = partMetric.cell_last_overlap_frame;
        mrow.cell_overlap_duration_frames = partMetric.cell_overlap_duration_frames;
        mrow.cell_peak_overlap_npix = partMetric.cell_peak_overlap_npix;
        mrow.cell_peak_overlap_frac_of_cell = partMetric.cell_peak_overlap_frac_of_cell;
        mrow.cell_peak_overlap_frac_of_event_frame = partMetric.cell_peak_overlap_frac_of_event_frame;
        if ~isnan(lmkIdx) && numel(lmkPixCache) >= lmkIdx
            lmkPix = lmkPixCache{lmkIdx};
            mrow.soma_npix = numel(lmkPix);
            if ~isempty(evtPix) && ~isempty(lmkPix)
                mrow.soma_overlap_npix = numel(intersect(evtPix, lmkPix));
                if mrow.soma_npix > 0
                    mrow.soma_overlap_frac_of_soma = mrow.soma_overlap_npix / mrow.soma_npix;
                end
                if mrow.event_spatial_npix > 0
                    mrow.soma_overlap_frac_of_event = mrow.soma_overlap_npix / mrow.event_spatial_npix;
                end
                mrow.soma_covered_ge50 = mrow.soma_overlap_frac_of_soma >= 0.5;
            end
            lmkCenter = [];
            if isnumeric(localLmkCenter) && size(localLmkCenter,1) >= lmkIdx
                lmkCenter = localLmkCenter(lmkIdx,:);
            end
            somaMask = [];
            if iscell(localLmkMask) && numel(localLmkMask) >= lmkIdx
                somaMask = localLmkMask{lmkIdx};
            end
            somaMetric = canonicalSomaDistanceMetrics(evtPix, xSpaTempCache{ee}, all_t0(idx), all_t1(idx), [sz_x sz_y sz_z sz_t], somaMask, lmkCenter);
            mrow.soma_centroid_dist_px = somaMetric.soma_centroid_dist_px;
            mrow.soma_start_dist_px = somaMetric.soma_start_dist_px;
            mrow.soma_end_dist_px = somaMetric.soma_end_dist_px;
            mrow.landmark_dist_avg = safe2d(all_lmkDistAvg, idx, lmkIdx);
            mrow.landmark_dist_min = safe2d(all_lmkDistMin, idx, lmkIdx);
            mrow.chg_toward = safe2d(all_chgToward, idx, lmkIdx);
            mrow.chg_away = safe2d(all_chgAway, idx, lmkIdx);
            mrow.chg_toward_bef_reach = safe2d(all_chgTowardBef, idx, lmkIdx);
            mrow.chg_away_aft_reach = safe2d(all_chgAwayAft, idx, lmkIdx);
            den = mrow.chg_toward + mrow.chg_away;
            if isfinite(den) && den > 0
                mrow.motion_index = (mrow.chg_toward - mrow.chg_away) / den;
            end
            mrow.motion_class = motionClass(mrow.chg_toward, mrow.chg_away, mrow.chg_toward_bef_reach, mrow.chg_away_aft_reach);
        else
            mrow.motion_class = "no_landmark";
        end
        eventMapRows(mm,1) = mrow;
    end
    primaryEntry = NaN;
    primaryEntryThresh = NaN;
    primaryMask = [eventMapRows.is_primary_cell];
    if any(primaryMask)
        primaryRow = eventMapRows(find(primaryMask, 1));
        primaryEntry = primaryRow.cell_first_overlap_frame;
        if isfield(primaryRow, 'cell_entry_frame_cov_thresh')
            primaryEntryThresh = primaryRow.cell_entry_frame_cov_thresh;
        end
    end
    if isfinite(primaryEntry)
        for mm = 1:numel(eventMapRows)
            if isfinite(eventMapRows(mm).cell_first_overlap_frame)
                eventMapRows(mm).cell_entry_delay_from_primary_frames = eventMapRows(mm).cell_first_overlap_frame - primaryEntry;
            end
        end
    end
    if ~isfinite(primaryEntryThresh)
        primaryEntryThresh = primaryEntry;
    end
    if isfinite(primaryEntryThresh)
        for mm = 1:numel(eventMapRows)
            cellEntryThresh = eventMapRows(mm).cell_entry_frame_cov_thresh;
            if isfinite(cellEntryThresh)
                eventMapRows(mm).cell_entry_delay_from_primary_cov_thresh_frames = cellEntryThresh - primaryEntryThresh;
            elseif eventMapRows(mm).is_primary_cell && isfinite(primaryEntryThresh)
                eventMapRows(mm).cell_entry_delay_from_primary_cov_thresh_frames = 0;
            end
        end
    end
    eventMapRows = localFilterOfficialCellMembers(eventMapRows, minAssociatedCellCoverageFrac);
    mapRows = [mapRows; eventMapRows]; %#ok<AGROW>
end
end

function out = fileTemplate()
out = struct('file_id',"",'animal_id',"",'slice_num',"",'filepath',"",'filename',"",'exp_date',"",...
    'sz_x',NaN,'sz_y',NaN,'sz_z',NaN,'sz_t',NaN,...
    'frame_rate_s',NaN,'n_frames',NaN,'spatial_res',NaN,'n_events',NaN,'n_cells',NaN,'n_somas',NaN,...
    'animal_name',"",'genotype',"",'age',"",'gender',"",'region_imaged',"",'virus',"",'drug',"",'stim',"",'training',"",'image_epoc',"");
end
function out = cellTemplate()
out = struct('cell_uid',"",'file_id',"",'animal_id',"",'slice_num',"",'local_cell_id',NaN,'is_soma_like',false,'soma_landmark_idx',NaN,'cell_mask_npix',NaN);
end
function out = eventTemplate()
out = struct('event_uid',"",'file_id',"",'animal_id',"",'slice_num',"",'local_event_id',NaN,'t0_frame',NaN,'t1_frame',NaN,'t0_sec',NaN,'t1_sec',NaN,...
    'area_px2',NaN,'duration_frames',NaN,'rise19_frames',NaN,'fall91_frames',NaN,'decay_tau_frames',NaN,'dff_auc',NaN,'dff_max2',NaN,'max_prop_speed',NaN,...
    'perimeter_px',NaN,'circularity',NaN,'propagation_distance_px',NaN,'propagation_angle_deg',NaN,...
    'branch_origin_x_px',NaN,'branch_origin_y_px',NaN,'branch_origin_frame',NaN,'branch_origin_cell_uid',"",...
    'primary_cell_uid',"");
end
function out = mapTemplate()
out = struct('event_uid',"",'cell_uid',"",'association_score',NaN,'is_primary_cell',false, ...
    'soma_landmark_idx',NaN,'event_spatial_npix',NaN,'cell_overlap_npix',NaN,'cell_overlap_frac_of_event',NaN,...
    'cell_overlap_frac_of_cell',NaN,'frac_outside_cell',NaN,'outside_cell_area_npix',NaN,'soma_npix',NaN,'soma_overlap_npix',NaN, ...
    'cell_first_overlap_frame',NaN,'cell_entry_frame_cov_thresh',NaN,'cell_last_overlap_frame',NaN,'cell_overlap_duration_frames',NaN, ...
    'cell_peak_overlap_npix',NaN,'cell_peak_overlap_frac_of_cell',NaN,'cell_peak_overlap_frac_of_event_frame',NaN, ...
    'cell_entry_delay_from_primary_frames',NaN,'cell_entry_delay_from_primary_cov_thresh_frames',NaN, ...
    'soma_overlap_frac_of_soma',NaN,'soma_overlap_frac_of_event',NaN,'soma_covered_ge50',false, ...
    'soma_centroid_dist_px',NaN,'soma_start_dist_px',NaN,'soma_end_dist_px',NaN, ...
    'landmark_dist_avg',NaN,'landmark_dist_min',NaN,'chg_toward',NaN,'chg_away',NaN, ...
    'chg_toward_bef_reach',NaN,'chg_away_aft_reach',NaN,'motion_index',NaN,'motion_class',"");
end

function meta = findFileMeta(DB_Files, fileID)
meta = struct('filepath','','filename','','expDate','');
if isempty(DB_Files), return; end
idx = find(strcmp({DB_Files.fileID}, fileID), 1);
if ~isempty(idx), meta = DB_Files(idx); end
end
function meta = findAnimalMeta(DB_Animals, animalID)
meta = struct('animalName','','genotype','','age','','gender','');
if isempty(DB_Animals) || isempty(animalID), return; end
idx = find(strcmp({DB_Animals.animalID}, animalID), 1);
if ~isempty(idx), meta = DB_Animals(idx); end
end
function meta = findExpMeta(DB_ExpCond, fileID)
meta = struct('regionImaged','','virus','','drug','','stim','','training','','imageEpoc','');
if isempty(DB_ExpCond) || isempty(fileID), return; end
idx = find(strcmp({DB_ExpCond.fileID}, fileID), 1);
if isempty(idx), return; end
row = DB_ExpCond(idx);
meta.regionImaged = getFieldOr(row, 'regionImaged', '');
meta.virus = getFieldOr(row, 'virus', '');
meta.drug = getFieldOr(row, 'drug', '');
meta.stim = getFieldOr(row, 'stim', '');
meta.training = getFieldOr(row, 'training', '');
meta.imageEpoc = getFieldOr(row, 'imageEpoc', '');
end

function out = getFieldOr(s, f, d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), out = s.(f); else, out = d; end
end
function out = getFieldPathOr(s, dotted, d)
out = d; cur = s; parts = strsplit(dotted,'.');
for ii=1:numel(parts)
    if ~isstruct(cur) || ~isfield(cur,parts{ii}), return; end
    cur = cur.(parts{ii});
end
out = cur;
end
function c = toCellstr(v)
if isempty(v), c = {}; return; end
if iscell(v), c = cellfun(@toStr,v,'UniformOutput',false); return; end
if isstring(v), c = cellstr(v(:)); return; end
if ischar(v), c = {v}; return; end
c = cellstr(string(v(:)));
end
function s = toStr(v)
if ischar(v), s=v; elseif isstring(v), s=char(v); elseif isempty(v), s=''; else, s=char(string(v)); end
end
function c = asChar(v), c = toStr(v); end
function v = toScalarOrNaN(x, idx)
v = NaN;
if isempty(x), return; end
if iscell(x)
    if numel(x)>=idx, xv = x{idx}; if isnumeric(xv)&&~isempty(xv), v=xv(1); end, end
else
    if numel(x)>=idx, xv = x(idx); if isnumeric(xv), v=xv; end, end
end
end
function v = toScalarOrEmpty(x, idx)
v = '';
if isempty(x), return; end
if iscell(x)
    if numel(x)>=idx, v = x{idx}; end
else
    if numel(x)>=idx, v = x(idx); end
end
end
function v = forceCol(x)
if isempty(x), v=[]; else, v=x(:); end
end
function v = safeIndex(x, idx)
v = NaN;
if isempty(x), return; end
if numel(x)>=idx, v = x(idx); end
end

function v = safeCellIndex(c, idx)
v = [];
if iscell(c) && numel(c) >= idx
    v = c{idx};
end
end

function [primaryCell, mappedCells, mappedScores] = extractMemberMapping(memberIdx, idx)
primaryCell = NaN;
mappedCells = [];
mappedScores = [];
if isempty(memberIdx), return; end

if isnumeric(memberIdx)
    if isvector(memberIdx)
        if numel(memberIdx) >= idx
            cid = memberIdx(idx);
            if ~isnan(cid) && cid > 0
                primaryCell = cid;
                mappedCells = cid;
                mappedScores = NaN;
            end
        end
        return
    end
    if size(memberIdx,1) >= idx
        row = memberIdx(idx,:);
        valid = find(row > 0);
        if ~isempty(valid)
            primaryCell = valid(1);
            mappedCells = valid(:)';
            mappedScores = row(valid);
        end
        return
    end
end

if iscell(memberIdx) && numel(memberIdx) >= idx
    ids = memberIdx{idx};
    if isempty(ids), return; end
    ids = ids(:);
    ids = ids(~isnan(ids) & ids > 0);
    if isempty(ids), return; end
    primaryCell = ids(1);
    mappedCells = ids(:)';
    mappedScores = NaN(size(mappedCells));
end
end

function [mappedCells, mappedScores] = localMergePrimaryIntoCandidates(mappedCells, mappedScores, primaryCell)
mappedCells = double(mappedCells(:)');
mappedCells = mappedCells(isfinite(mappedCells) & mappedCells > 0);
mappedScores = double(mappedScores(:)');
if numel(mappedScores) < numel(mappedCells)
    mappedScores(numel(mappedScores)+1:numel(mappedCells)) = NaN;
elseif numel(mappedScores) > numel(mappedCells)
    mappedScores = mappedScores(1:numel(mappedCells));
end
if ~(isfinite(primaryCell) && primaryCell > 0) || any(mappedCells == primaryCell)
    return
end
mappedCells = [primaryCell, mappedCells];
mappedScores = [NaN, mappedScores];
end

function rows = localFilterOfficialCellMembers(rows, minCoverageFrac)
if isempty(rows)
    return
end
if ~isfinite(minCoverageFrac)
    minCoverageFrac = 0.10;
end
isPrimary = [rows.is_primary_cell]';
peakFrac = NaN(numel(rows), 1);
for ii = 1:numel(rows)
    if isfield(rows(ii), 'cell_peak_overlap_frac_of_cell')
        peakFrac(ii) = double(rows(ii).cell_peak_overlap_frac_of_cell);
    end
end
keep = isPrimary | (isfinite(peakFrac) & peakFrac >= minCoverageFrac);
rows = rows(keep);
end

function ownerCell = localResolveOwnerCellFromPoint(pt, cellMasks, candidateCells)
ownerCell = NaN;
if any(~isfinite(pt)) || ~iscell(cellMasks) || isempty(cellMasks)
    return
end
row = round(pt(2));
col = round(pt(1));
if isempty(candidateCells)
    candidateCells = 1:numel(cellMasks);
end
candidateCells = double(candidateCells(:)');
candidateCells = candidateCells(isfinite(candidateCells) & candidateCells >= 1 & candidateCells <= numel(cellMasks));
inside = [];
for ii = 1:numel(candidateCells)
    cc = candidateCells(ii);
    mask = logical(cellMasks{cc});
    if isempty(mask) || row < 1 || col < 1 || row > size(mask,1) || col > size(mask,2)
        continue
    end
    if mask(row, col)
        inside(end+1) = cc; %#ok<AGROW>
    end
end
if isempty(inside)
    ownerCell = localNearestCellMaskCentroid(pt, candidateCells, cellMasks);
else
    ownerCell = localNearestCellMaskCentroid(pt, inside, cellMasks);
end
end

function ownerCell = localNearestCellMaskCentroid(pt, candidateCells, cellMasks)
ownerCell = NaN;
bestD = inf;
for ii = 1:numel(candidateCells)
    cc = candidateCells(ii);
    mask = logical(cellMasks{cc});
    [r, c] = find(mask);
    if isempty(r)
        continue
    end
    ctr = [mean(c), mean(r)];
    d = hypot(ctr(1) - pt(1), ctr(2) - pt(2));
    if d < bestD
        bestD = d;
        ownerCell = cc;
    end
end
end

function p = getEventSpatialPixels(xSpa, xSpaTemp, idx, opts, dbIdx)
p = [];
if iscell(xSpa) && numel(xSpa) >= idx && isnumeric(xSpa{idx})
    p = unique(xSpa{idx}(:));
    p = p(p > 0 & isfinite(p));
    if ~isempty(p)
        return
    end
end
if iscell(xSpaTemp) && numel(xSpaTemp) >= idx && isnumeric(xSpaTemp{idx})
    p = voxelToSpatial(xSpaTemp{idx}, opts, dbIdx);
end
end

function p = voxelToSpatial(vIdx, opts, dbIdx)
p = [];
if ~isfield(opts, 'sz') || size(opts.sz,1) < dbIdx || size(opts.sz,2) < 4
    return
end
sx = opts.sz(dbIdx,1); sy = opts.sz(dbIdx,2); sz = opts.sz(dbIdx,3); st = opts.sz(dbIdx,4);
if any(~isfinite([sx sy sz st])) || any([sx sy sz st] <= 0)
    return
end
try
    [ih, iw, il, ~] = ind2sub([sx sy sz st], vIdx(:));
    p = unique(sub2ind([sx sy sz], ih, iw, il));
catch
    p = [];
end
end

function p = maskToPix(m)
p = [];
if isnumeric(m) || islogical(m)
    p = find(m > 0);
end
end

function cache = cacheMaskPixels(masks)
cache = {};
if ~iscell(masks)
    return
end
cache = cell(numel(masks), 1);
for ii = 1:numel(masks)
    cache{ii} = maskToPix(masks{ii});
end
end

function out = sliceCellBlock(cellsIn, offset, nLocal)
out = {};
if ~iscell(cellsIn) || nLocal <= 0
    return
end
idx1 = offset + 1;
idx2 = min(numel(cellsIn), offset + nLocal);
if idx1 > idx2
    return
end
out = cellsIn(idx1:idx2);
end

function out = sliceNumericBlock(x, offset, nRowsLocal, nColsLocal)
out = [];
if ~isnumeric(x) || isempty(x) || nRowsLocal <= 0 || nColsLocal <= 0
    return
end
row1 = offset + 1;
row2 = min(size(x,1), offset + nRowsLocal);
col2 = min(size(x,2), nColsLocal);
if row1 > row2 || col2 < 1
    return
end
out = x(row1:row2, 1:col2);
end

function out = sliceNumericRows(x, offset, nRowsLocal)
out = [];
if ~isnumeric(x) || isempty(x) || nRowsLocal <= 0
    return
end
row1 = offset + 1;
row2 = min(size(x,1), offset + nRowsLocal);
if row1 > row2
    return
end
out = x(row1:row2, :);
end

function n = scalarCount(v)
n = 0;
if isnumeric(v) && isfinite(v) && v > 0
    n = round(v);
end
end

function v = safe2d(x, i, j)
v = NaN;
if isempty(x) || ~isnumeric(x)
    return
end
if size(x,1) >= i && size(x,2) >= j
    v = x(i,j);
end
end

function cls = motionClass(toward, away, towardBeforeReach, awayAfterReach)
cls = "undetermined";
if isfinite(towardBeforeReach) && towardBeforeReach > 0 && isfinite(awayAfterReach) && awayAfterReach > 0
    cls = "traverse";
    return
end
if isfinite(toward) && isfinite(away)
    if toward > away
        cls = "toward";
    elseif away > toward
        cls = "away";
    elseif toward == 0 && away == 0
        cls = "static";
    else
        cls = "mixed";
    end
end
end

function cellToSoma = resolveSingleSomaPerCell(inclu, cellMasks, lmkMasks, cellCenters, lmkCenters)
if ~isnumeric(inclu) || isempty(inclu)
    cellToSoma = NaN(0,1);
    return
end
nCells = size(inclu,1);
cellToSoma = NaN(nCells,1);
for cc = 1:nCells
    cand = find(inclu(cc,:) > 0);
    if isempty(cand)
        continue
    end
    if numel(cand) == 1
        cellToSoma(cc) = cand(1);
        continue
    end
    cellToSoma(cc) = pickMostCentralLandmark(cand, cc, cellMasks, lmkMasks, cellCenters, lmkCenters);
end
end

function best = pickMostCentralLandmark(cand, cellIdx, cellMasks, lmkMasks, cellCenters, lmkCenters)
best = cand(1);
scoreOverlap = -inf(numel(cand),1);
scoreDist = inf(numel(cand),1);
cMask = {};
if iscell(cellMasks) && numel(cellMasks) >= cellIdx
    cMask = cellMasks(cellIdx);
end
cCenter = getCenterRow(cellCenters, cellIdx);
if any(isnan(cCenter)) && ~isempty(cMask)
    cCenter = maskCenter(cMask{1});
end
for ii = 1:numel(cand)
    ll = cand(ii);
    lm = [];
    if iscell(lmkMasks) && numel(lmkMasks) >= ll
        lm = lmkMasks{ll};
    end
    if ~isempty(cMask) && ~isempty(lm) && (isnumeric(lm) || islogical(lm))
        cm = cMask{1};
        if isnumeric(cm) || islogical(cm)
            scoreOverlap(ii) = nnz((cm > 0) & (lm > 0));
        end
    end
    lCenter = getCenterRow(lmkCenters, ll);
    if any(isnan(lCenter)) && ~isempty(lm)
        lCenter = maskCenter(lm);
    end
    if all(isfinite(cCenter(1:2))) && all(isfinite(lCenter(1:2)))
        scoreDist(ii) = hypot(cCenter(1)-lCenter(1), cCenter(2)-lCenter(2));
    end
end
[~, ord] = sortrows([ -scoreOverlap, scoreDist, cand(:) ], [1 2 3]);
best = cand(ord(1));
end

function c = getCenterRow(C, idx)
c = [NaN NaN NaN];
if isnumeric(C) && size(C,1) >= idx
    row = C(idx,:);
    c(1:min(3,numel(row))) = row(1:min(3,numel(row)));
end
end

function c = maskCenter(m)
c = [NaN NaN NaN];
if isempty(m) || ~(isnumeric(m) || islogical(m))
    return
end
nd = ndims(m);
if nd <= 2
    [r, col] = find(m > 0);
    if ~isempty(r)
        c = [mean(r), mean(col), 1];
    end
else
    [r, col, z] = ind2sub(size(m), find(m > 0));
    if ~isempty(r)
        c = [mean(r), mean(col), mean(z)];
    end
end
end
