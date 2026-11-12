function out = ingestAQuA2Canonical(varargin)
    % Canonical ingest for AQUA2 outputs:
    % writes immutable run-stamped tables: files, cells, events, event_cell_map.

    BWA_setupPaths(fileparts(mfilename('fullpath')));

    p = inputParser;
    p.StructExpand = false;
    addParameter(p, 'rootPath', '.', @ischar);
    addParameter(p, 'outputRoot', fullfile('derived', 'ingest_runs'), @ischar);
    addParameter(p, 'saveCsv', false, @islogical);
    addParameter(p, 'useParallel', false, @islogical);
    addParameter(p, 'minAssociatedCellCoverageFrac', 0.10, @isnumeric);
    parse(p, varargin{:});

    rootPath = p.Results.rootPath;
    outputRoot = p.Results.outputRoot;
    saveCsv = p.Results.saveCsv;
    useParallel = p.Results.useParallel;
    minAssociatedCellCoverageFrac = double(p.Results.minAssociatedCellCoverageFrac);

    tAll = tic;
    disp('ingestAQuA2Canonical: starting');
    disp(['  rootPath=' rootPath ', outputRoot=' outputRoot ', saveCsv=' tf(saveCsv) ', useParallel=' tf(useParallel)]);

    savedPwd = pwd;
    c = onCleanup(@() cd(savedPwd));
    cd(rootPath);

    disp('Canonical raw: building databases...');
    [DB_Files, DB_Animals, DB_ExpCond] = CreateTopDataBases('saveLegacyDbs', false);

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    runDir = fullfile(outputRoot, timestamp);
    if ~exist(runDir, 'dir')
        mkdir(runDir);
    end

    nF = numel(DB_Files);
    disp(['Canonical raw: processing ' num2str(nF) ' files']);
    fileChunks = cell(nF,1);
    cellChunks = cell(nF,1);
    eventChunks = cell(nF,1);
    mapChunks = cell(nF,1);
    skipPaths = cell(nF,1);
    skipReasons = cell(nF,1);

    didParallel = false;
    if useParallel
        try
            if isempty(gcp('nocreate'))
                parpool;
            end
            pObj = gcp('nocreate');
            disp(['Canonical raw: using parallel pool with ' num2str(pObj.NumWorkers) ' workers']);
            q = parallel.pool.DataQueue;
            nDone = 0;
            reportStep = max(1, floor(nF/20));
            afterEach(q, @parallelTick);
            parfor ii = 1:nF
                [fileChunks{ii}, cellChunks{ii}, eventChunks{ii}, mapChunks{ii}, skipPaths{ii}, skipReasons{ii}] = ...
                    processOneFile(DB_Files(ii), DB_Animals, DB_ExpCond, minAssociatedCellCoverageFrac);
                send(q, ii);
            end
            didParallel = true;
        catch
            didParallel = false;
        end
    end

    if ~didParallel
        if useParallel
            disp('Canonical raw: parallel unavailable, falling back to sequential mode');
        else
            disp('Canonical raw: using sequential mode');
        end
        for ii = 1:nF
            disp(['Canonical seq [' num2str(ii) '/' num2str(nF) ']']);
            [fileChunks{ii}, cellChunks{ii}, eventChunks{ii}, mapChunks{ii}, skipPaths{ii}, skipReasons{ii}] = ...
                processOneFile(DB_Files(ii), DB_Animals, DB_ExpCond, minAssociatedCellCoverageFrac);
        end
    end

    disp('Canonical raw: merging chunks and writing outputs...');
    fileRows = repmat(makeFileRowTemplate(), 0, 1);
    cellRows = repmat(makeCellRowTemplate(), 0, 1);
    eventRows = repmat(makeEventRowTemplate(), 0, 1);
    mapRows = repmat(makeMapRowTemplate(), 0, 1);
    skipped = {};
    skippedReasons = {};
    for ii = 1:nF
        if ~isempty(fileChunks{ii}), fileRows = [fileRows; fileChunks{ii}]; end %#ok<AGROW>
        if ~isempty(cellChunks{ii}), cellRows = [cellRows; cellChunks{ii}]; end %#ok<AGROW>
        if ~isempty(eventChunks{ii}), eventRows = [eventRows; eventChunks{ii}]; end %#ok<AGROW>
        if ~isempty(mapChunks{ii}), mapRows = [mapRows; mapChunks{ii}]; end %#ok<AGROW>
        if ~isempty(skipPaths{ii})
            skipped{end+1,1} = skipPaths{ii}; %#ok<AGROW>
            skippedReasons{end+1,1} = skipReasons{ii}; %#ok<AGROW>
        end
    end

    files = struct2table(fileRows);
    cells = struct2table(cellRows);
    events = struct2table(eventRows);
    event_cell_map = struct2table(mapRows);

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
    manifest.root_path = asString(rootPath);
    manifest.run_dir = asString(runDir);
    manifest.n_files = height(files);
    manifest.n_cells = height(cells);
    manifest.n_events = height(events);
    manifest.n_event_cell_map = height(event_cell_map);
    manifest.skipped_files = skipped;
    manifest.skipped_reasons = skippedReasons;
    manifest.created_on = char(datetime('now'));

    save(fullfile(runDir, 'manifest.mat'), 'manifest');
    fid = fopen(fullfile(runDir, 'manifest.json'), 'w');
    if fid > 0
        fwrite(fid, jsonencode(manifest), 'char');
        fclose(fid);
    end

    latestPointer = fullfile(outputRoot, 'latest_ingest.txt');
    fid = fopen(latestPointer, 'w');
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

    disp(['ingestAQuA2Canonical: done in ' num2str(toc(tAll), '%.1f') ' s']);

    function parallelTick(~)
        nDone = nDone + 1;
        if nDone == 1 || nDone == nF || mod(nDone, reportStep) == 0
            disp(['Canonical par [' num2str(nDone) '/' num2str(nF) ']']);
        end
    end
end

function [fRows, cRows, eRows, mRows, skipPath, skipReason] = processOneFile(thisFile, DB_Animals, DB_ExpCond, minAssociatedCellCoverageFrac)
    fRows = repmat(makeFileRowTemplate(), 0, 1);
    cRows = repmat(makeCellRowTemplate(), 0, 1);
    eRows = repmat(makeEventRowTemplate(), 0, 1);
    mRows = repmat(makeMapRowTemplate(), 0, 1);
    skipPath = '';
    skipReason = '';

    fullPath = fullfile(thisFile.filepath, thisFile.filename);
    try
        mObj = matfile(fullPath);
        if ~ismember('res', who(mObj))
            skipPath = fullPath;
            skipReason = 'missing variable: res';
            return
        end
        resStruct = mObj.res;
        if ~isstruct(resStruct) || ~isfield(resStruct, 'fts1') || ~isfield(resStruct, 'opts')
            skipPath = fullPath;
            skipReason = 'res missing fts1/opts';
            return
        end
        fts1 = resStruct.fts1;
        opts = resStruct.opts;
    catch ME
        skipPath = fullPath;
        skipReason = ME.message;
        return
    end

    idMini = getIDfileBy(thisFile.filepath, 'filename', thisFile.filename);
    if isempty(idMini)
        idMini = struct('filename', thisFile.filename, 'fileID', thisFile.fileID, 'animalID', '', 'sliceNum', '');
    end

    animalMeta = getAnimalMeta(DB_Animals, idMini.animalID);
    expMeta = getExpMeta(DB_ExpCond, idMini.fileID);

    nEvents = numel(getPathOrDefault(fts1, 'loc.t0', []));
    nCells = inferNCells(fts1);
    nSomas = inferNSomas(fts1);

    fr = scalarOrNaN(getPathOrDefault(opts, 'frameRate', NaN));
    sz = getPathOrDefault(opts, 'sz', []);
    nFrames = NaN;
    sz_x = NaN; sz_y = NaN; sz_z = NaN; sz_t = NaN;
    if isnumeric(sz) && ~isempty(sz)
        if numel(sz) >= 1, sz_x = sz(1); end
        if numel(sz) >= 2, sz_y = sz(2); end
        if numel(sz) >= 3, sz_z = sz(3); end
        if numel(sz) >= 4, sz_t = sz(4); end
        if size(sz, 2) >= 4
            nFrames = sz(1, 4);
        elseif numel(sz) >= 4
            nFrames = sz(4);
        end
    end
    spatialRes = scalarOrNaN(getPathOrDefault(opts, 'spatialRes', NaN));

    fRow = makeFileRowTemplate();
    fRow.file_id = asString(idMini.fileID);
    fRow.animal_id = asString(idMini.animalID);
    fRow.slice_num = asString(idMini.sliceNum);
    fRow.filepath = asString(thisFile.filepath);
    fRow.filename = asString(thisFile.filename);
    fRow.exp_date = asString(thisFile.expDate);
    fRow.sz_x = sz_x;
    fRow.sz_y = sz_y;
    fRow.sz_z = sz_z;
    fRow.sz_t = sz_t;
    fRow.frame_rate_s = fr;
    fRow.n_frames = nFrames;
    fRow.spatial_res = spatialRes;
    fRow.n_events = nEvents;
    fRow.n_cells = nCells;
    fRow.n_somas = nSomas;
    fRow.animal_name = asString(animalMeta.animalName);
    fRow.genotype = asString(animalMeta.genotype);
    fRow.age = asString(animalMeta.age);
    fRow.gender = asString(animalMeta.gender);
    fRow.region_imaged = asString(expMeta.regionImaged);
    fRow.virus = asString(expMeta.virus);
    fRow.drug = asString(expMeta.drug);
    fRow.stim = asString(expMeta.stim);
    fRow.training = asString(expMeta.training);
    fRow.image_epoc = asString(expMeta.imageEpoc);
    fRows = fRow;

    fileID = asChar(idMini.fileID);
    animalID = asChar(idMini.animalID);
    sliceNum = asChar(idMini.sliceNum);
    cRows = buildCellRows(fts1, fileID, animalID, sliceNum, nCells);
    [eRows, mRows] = buildEventRows(fts1, fileID, animalID, sliceNum, fr, spatialRes, [sz_x sz_y sz_z sz_t], minAssociatedCellCoverageFrac);
end

function out = makeFileRowTemplate()
    out = struct( ...
        'file_id', "", 'animal_id', "", 'slice_num', "", ...
        'filepath', "", 'filename', "", 'exp_date', "", ...
        'sz_x', NaN, 'sz_y', NaN, 'sz_z', NaN, 'sz_t', NaN, ...
        'frame_rate_s', NaN, 'n_frames', NaN, 'spatial_res', NaN, ...
        'n_events', NaN, 'n_cells', NaN, 'n_somas', NaN, ...
        'animal_name', "", 'genotype', "", 'age', "", 'gender', "", ...
        'region_imaged', "", 'virus', "", 'drug', "", 'stim', "", 'training', "", ...
        'image_epoc', "");
end

function out = makeCellRowTemplate()
    out = struct( ...
        'cell_uid', "", 'file_id', "", 'animal_id', "", 'slice_num', "", ...
        'local_cell_id', NaN, 'is_soma_like', false, 'soma_landmark_idx', NaN, 'cell_mask_npix', NaN);
end

function out = makeEventRowTemplate()
    out = struct( ...
        'event_uid', "", 'file_id', "", 'animal_id', "", 'slice_num', "", ...
        'local_event_id', NaN, 't0_frame', NaN, 't1_frame', NaN, 't0_sec', NaN, 't1_sec', NaN, ...
        'area_px2', NaN, 'duration_frames', NaN, 'rise19_frames', NaN, 'fall91_frames', NaN, 'decay_tau_frames', NaN, ...
        'dff_auc', NaN, 'dff_max2', NaN, 'max_prop_speed', NaN, ...
        'perimeter_px', NaN, 'circularity', NaN, 'propagation_distance_px', NaN, ...
        'propagation_angle_deg', NaN, ...
        'branch_origin_x_px', NaN, 'branch_origin_y_px', NaN, 'branch_origin_frame', NaN, 'branch_origin_cell_uid', "", ...
        'primary_cell_uid', "");
end

function out = makeMapRowTemplate()
    out = struct( ...
        'event_uid', "", 'cell_uid', "", ...
        'association_score', NaN, 'is_primary_cell', false, ...
        'soma_landmark_idx', NaN, 'event_spatial_npix', NaN, ...
        'cell_overlap_npix', NaN, 'cell_overlap_frac_of_event', NaN, ...
        'cell_overlap_frac_of_cell', NaN, 'frac_outside_cell', NaN, 'outside_cell_area_npix', NaN, ...
        'cell_first_overlap_frame', NaN, 'cell_entry_frame_cov_thresh', NaN, 'cell_last_overlap_frame', NaN, 'cell_overlap_duration_frames', NaN, ...
        'cell_peak_overlap_npix', NaN, 'cell_peak_overlap_frac_of_cell', NaN, 'cell_peak_overlap_frac_of_event_frame', NaN, ...
        'cell_entry_delay_from_primary_frames', NaN, 'cell_entry_delay_from_primary_cov_thresh_frames', NaN, ...
        'soma_npix', NaN, 'soma_overlap_npix', NaN, ...
        'soma_overlap_frac_of_soma', NaN, 'soma_overlap_frac_of_event', NaN, 'soma_covered_ge50', false, ...
        'soma_centroid_dist_px', NaN, 'soma_start_dist_px', NaN, 'soma_end_dist_px', NaN, ...
        'landmark_dist_avg', NaN, 'landmark_dist_min', NaN, ...
        'chg_toward', NaN, 'chg_away', NaN, 'chg_toward_bef_reach', NaN, 'chg_away_aft_reach', NaN, ...
        'motion_index', NaN, 'motion_class', "");
end

function DB_ExpCond = loadExpCondIfPresent()
    if exist('DBExpCond.mat', 'file')
        loaded = load('DBExpCond.mat');
        if isfield(loaded, 'DB_ExpCond')
            DB_ExpCond = loaded.DB_ExpCond;
            return
        end
    end
    DB_ExpCond = struct('fileID',{},'animalID',{},'regionImaged',{},...
        'virus',{},'virusLoc',{},'dpi_virus',{},...
        'drug',{}, 'drugLoc', {}, 'dpi_drug',{},...
        'stim', {}, 'stimloc', {}, 'dpi_stim', {},...
        'training', {}, 'dpi_training', {},...
        'imageEpoc', {}, 'epocFrameStart', {}, 'epocFrameEnd', {}, 'epocTimeStart', {}, 'epocTimeEnd',{});
end

function meta = getAnimalMeta(DB_Animals, animalID)
    meta = struct('animalName', '', 'genotype', '', 'age', '', 'gender', '');
    if isempty(DB_Animals) || isempty(animalID)
        return
    end
    idx = find(strcmp({DB_Animals.animalID}, asChar(animalID)), 1);
    if isempty(idx)
        return
    end
    meta = DB_Animals(idx);
end

function meta = getExpMeta(DB_ExpCond, fileID)
    meta = struct('regionImaged', '', 'virus', '', 'drug', '', 'stim', '', 'training', '', 'imageEpoc', '');
    if isempty(DB_ExpCond) || isempty(fileID)
        return
    end
    idx = find(strcmp({DB_ExpCond.fileID}, asChar(fileID)), 1);
    if isempty(idx)
        return
    end
    row = DB_ExpCond(idx);
    meta.regionImaged = pickFieldOrEmpty(row, 'regionImaged');
    meta.virus = pickFieldOrEmpty(row, 'virus');
    meta.drug = pickFieldOrEmpty(row, 'drug');
    meta.stim = pickFieldOrEmpty(row, 'stim');
    meta.training = pickFieldOrEmpty(row, 'training');
    meta.imageEpoc = pickFieldOrEmpty(row, 'imageEpoc');
end

function out = pickFieldOrEmpty(s, fname)
    out = '';
    if isfield(s, fname) && ~isempty(s.(fname))
        out = s.(fname);
    end
end

function nCells = inferNCells(fts1)
    nCells = NaN;
    inclu = getPathOrDefault(fts1, 'region.cell.incluLmk', []);
    if isnumeric(inclu) && ~isempty(inclu)
        nCells = size(inclu, 1);
        return
    end
    memberIdx = getPathOrDefault(fts1, 'region.cell.memberIdx', []);
    if isnumeric(memberIdx) && ~isempty(memberIdx)
        nCells = max(memberIdx(:));
        return
    end
    masks = getPathOrDefault(fts1, 'region.cell.mask', []);
    if iscell(masks)
        nCells = numel(masks);
    end
end

function nSomas = inferNSomas(fts1)
    nSomas = NaN;
    inclu = getPathOrDefault(fts1, 'region.cell.incluLmk', []);
    if isnumeric(inclu) && ~isempty(inclu)
        nSomas = size(inclu, 2);
    end
end

function cRows = buildCellRows(fts1, fileID, animalID, sliceNum, nCells)
    cRows = repmat(makeCellRowTemplate(), 0, 1);
    if ~(isnumeric(nCells) && ~isnan(nCells) && nCells > 0)
        return
    end

    masks = getPathOrDefault(fts1, 'region.cell.mask', []);
    inclu = getPathOrDefault(fts1, 'region.cell.incluLmk', []);
    lmkMasks = getPathOrDefault(fts1, 'region.landMark.mask', {});
    cellCenters = getPathOrDefault(fts1, 'region.cell.center', []);
    lmkCenters = getPathOrDefault(fts1, 'region.landMark.center', []);
    cellToSoma = resolveSingleSomaPerCell(inclu, masks, lmkMasks, cellCenters, lmkCenters);
    for cc = 1:nCells
        row = makeCellRowTemplate();
        row.file_id = asString(fileID);
        row.animal_id = asString(animalID);
        row.slice_num = asString(sliceNum);
        row.local_cell_id = cc;
        row.cell_uid = string(sprintf('%s_C%d', fileID, cc));
        row.is_soma_like = false;
        row.soma_landmark_idx = NaN;
        if cc <= numel(cellToSoma)
            row.soma_landmark_idx = cellToSoma(cc);
            row.is_soma_like = isfinite(row.soma_landmark_idx) && row.soma_landmark_idx > 0;
        end
        if iscell(masks) && numel(masks) >= cc
            thisMask = masks{cc};
            if isnumeric(thisMask) || islogical(thisMask)
                row.cell_mask_npix = nnz(thisMask);
            end
        end
        cRows(end+1,1) = row; %#ok<AGROW>
    end
end

function [eRows, mRows] = buildEventRows(fts1, fileID, animalID, sliceNum, frameRate, spatialRes, spatialSize, minAssociatedCellCoverageFrac)
    eRows = repmat(makeEventRowTemplate(), 0, 1);
    mRows = repmat(makeMapRowTemplate(), 0, 1);

    t0 = forceColumn(getPathOrDefault(fts1, 'loc.t0', []));
    if isempty(t0)
        return
    end
    nEvents = numel(t0);

    t1 = forceColumnWithLength(getPathOrDefault(fts1, 'loc.t1', []), nEvents);
    area = forceColumnWithLength(getPathOrDefault(fts1, 'basic.area', []), nEvents);
    duration = forceColumnWithLength(getPathOrDefault(fts1, 'curve.duration', []), nEvents);
    rise19 = forceColumnWithLength(getPathOrDefault(fts1, 'curve.rise19', []), nEvents);
    fall91 = forceColumnWithLength(getPathOrDefault(fts1, 'curve.fall91', []), nEvents);
    decayTau = forceColumnWithLength(getPathOrDefault(fts1, 'curve.decayTau', []), nEvents);
    dffAUC = forceColumnWithLength(getPathOrDefault(fts1, 'curve.dffAUC', []), nEvents);
    dffMax2 = forceColumnWithLength(getPathOrDefault(fts1, 'curve.dffMax2', []), nEvents);
    maxPropSpeed = forceColumnWithLength(getPathOrDefault(fts1, 'propagation.maxPropSpeed', []), nEvents);
    propGrow = getPathOrDefault(fts1, 'propagation.propGrow', {});
    propShrink = getPathOrDefault(fts1, 'propagation.propShrink', {});

    [mapPairs, primaryCellID] = eventCellPairs(fts1, nEvents);
    somaCtx = buildSomaContext(fts1, nEvents, spatialSize);
    somaCoverThr = 0.5;

    branchOriginXY = NaN(nEvents, 2);
    branchOriginFrame = NaN(nEvents, 1);
    branchOriginCellID = NaN(nEvents, 1);
    for ee = 1:nEvents
        theseCells = mapPairs(mapPairs(:,1) == ee, 2);
        tmpBWA = struct('fts1', fts1, 'opts', struct('sz', reshape(double(somaCtx.spatialSize), 1, []), ...
            'spatialRes', double(spatialRes), 'frameRate', double(frameRate)));
        owner = canonicalBranchOriginOwnership(tmpBWA, ee, ...
            'candidateCells', theseCells, ...
            'cellMasks', somaCtx.cellMasks);
        ownerCell = owner.ownerCell;
        originXY = owner.originXY;
        originFrame = owner.originFrame;
        if isfinite(ownerCell) && ownerCell > 0
            primaryCellID(ee) = ownerCell;
            branchOriginCellID(ee) = ownerCell;
            branchOriginXY(ee, :) = originXY;
            branchOriginFrame(ee) = originFrame;
        else
            origin = canonicalOriginCell(safeCellFetch(somaCtx.xSpaTemp, ee), t0(ee), t1(ee), ...
                somaCtx.spatialSize, theseCells, somaCtx.cellMasks);
            if isfinite(origin.originCell) && origin.originCell > 0
                primaryCellID(ee) = origin.originCell;
            end
        end
    end

    for ee = 1:nEvents
        row = makeEventRowTemplate();
        row.file_id = asString(fileID);
        row.animal_id = asString(animalID);
        row.slice_num = asString(sliceNum);
        row.local_event_id = ee;
        row.event_uid = string(sprintf('%s_E%d', fileID, ee));
        row.t0_frame = t0(ee);
        row.t1_frame = t1(ee);
        row.t0_sec = t0(ee) * frameRate;
        row.t1_sec = t1(ee) * frameRate;
        row.area_px2 = area(ee);
        row.duration_frames = duration(ee);
        row.rise19_frames = rise19(ee);
        row.fall91_frames = fall91(ee);
        row.decay_tau_frames = decayTau(ee);
        row.dff_auc = dffAUC(ee);
        row.dff_max2 = dffMax2(ee);
        row.max_prop_speed = maxPropSpeed(ee);
        evtPix = eventSpatialPixels(somaCtx, ee);
        evMetric = canonicalEventMetrics(evtPix, somaCtx.spatialSize, [], ...
            safeCellFetch(propGrow, ee), safeCellFetch(propShrink, ee));
        row.perimeter_px = evMetric.perimeter_px;
        row.circularity = evMetric.circularity;
        row.propagation_distance_px = evMetric.propagation_distance_px;
        row.propagation_angle_deg = evMetric.propagation_angle_deg;
        row.branch_origin_x_px = branchOriginXY(ee, 1);
        row.branch_origin_y_px = branchOriginXY(ee, 2);
        row.branch_origin_frame = branchOriginFrame(ee);
        if isfinite(branchOriginCellID(ee)) && branchOriginCellID(ee) > 0
            row.branch_origin_cell_uid = string(sprintf('%s_C%d', fileID, branchOriginCellID(ee)));
        else
            row.branch_origin_cell_uid = "";
        end
        if ~isnan(primaryCellID(ee)) && primaryCellID(ee) > 0
            row.primary_cell_uid = string(sprintf('%s_C%d', fileID, primaryCellID(ee)));
        else
            row.primary_cell_uid = "";
        end
        eRows(end+1,1) = row; %#ok<AGROW>
    end

    mapPairs = localEnsurePrimaryPairs(mapPairs, primaryCellID);
    eventMapRows = repmat(makeMapRowTemplate(), size(mapPairs, 1), 1);
    partitionFrameCache = cell(nEvents, 1);
    for ee = 1:nEvents
        partitionFrameCache{ee} = canonicalEventFrameCache(safeCellFetch(somaCtx.xSpaTemp, ee), ...
            t0(ee), t1(ee), somaCtx.spatialSize);
    end
    for rr = 1:size(mapPairs, 1)
        ev = mapPairs(rr,1);
        cc = mapPairs(rr,2);
        sc = mapPairs(rr,3);
        mp = makeMapRowTemplate();
        mp.event_uid = string(sprintf('%s_E%d', fileID, ev));
        mp.cell_uid = string(sprintf('%s_C%d', fileID, cc));
        mp.association_score = sc;
        mp.is_primary_cell = (~isnan(primaryCellID(ev)) && cc == primaryCellID(ev));
        if cc > 0 && cc <= numel(somaCtx.cellToSoma)
            lmk = somaCtx.cellToSoma(cc);
        else
            lmk = NaN;
        end
        mp.soma_landmark_idx = lmk;
        evtPix = eventSpatialPixels(somaCtx, ev);
        mp.event_spatial_npix = numel(evtPix);
        cellPix = [];
        if isfield(somaCtx, 'cellPix') && numel(somaCtx.cellPix) >= cc
            cellPix = somaCtx.cellPix{cc};
        end
        mapMetric = canonicalEventMetrics(evtPix, somaCtx.spatialSize, cellPix, [], []);
        mp.cell_overlap_npix = mapMetric.cell_overlap_npix;
        mp.cell_overlap_frac_of_event = mapMetric.cell_overlap_frac_of_event;
        mp.cell_overlap_frac_of_cell = mapMetric.cell_overlap_frac_of_cell;
        mp.frac_outside_cell = mapMetric.frac_outside_cell;
        mp.outside_cell_area_npix = mapMetric.outside_cell_area_npix;
        partMetric = canonicalEventCellPartitionMetrics(partitionFrameCache{ev}, ...
            t0(ev), t1(ev), somaCtx.spatialSize, cellPix, minAssociatedCellCoverageFrac);
        mp.cell_first_overlap_frame = partMetric.cell_first_overlap_frame;
        mp.cell_entry_frame_cov_thresh = partMetric.cell_entry_frame_cov_thresh;
        mp.cell_last_overlap_frame = partMetric.cell_last_overlap_frame;
        mp.cell_overlap_duration_frames = partMetric.cell_overlap_duration_frames;
        mp.cell_peak_overlap_npix = partMetric.cell_peak_overlap_npix;
        mp.cell_peak_overlap_frac_of_cell = partMetric.cell_peak_overlap_frac_of_cell;
        mp.cell_peak_overlap_frac_of_event_frame = partMetric.cell_peak_overlap_frac_of_event_frame;
        if ~isnan(lmk) && lmk > 0 && lmk <= numel(somaCtx.lmkPix)
            lmkPix = somaCtx.lmkPix{lmk};
            mp.soma_npix = numel(lmkPix);
            if ~isempty(evtPix) && ~isempty(lmkPix)
                mp.soma_overlap_npix = numel(intersect(evtPix, lmkPix));
                if mp.soma_npix > 0
                    mp.soma_overlap_frac_of_soma = mp.soma_overlap_npix / mp.soma_npix;
                end
                if mp.event_spatial_npix > 0
                    mp.soma_overlap_frac_of_event = mp.soma_overlap_npix / mp.event_spatial_npix;
                end
                mp.soma_covered_ge50 = mp.soma_overlap_frac_of_soma >= somaCoverThr;
            end
            lmkCenter = [];
            if isfield(somaCtx, 'lmkCenters') && isnumeric(somaCtx.lmkCenters) && size(somaCtx.lmkCenters,1) >= lmk
                lmkCenter = somaCtx.lmkCenters(lmk,:);
            end
            somaMetric = canonicalSomaDistanceMetrics(evtPix, safeCellFetch(somaCtx.xSpaTemp, ev), ...
                t0(ev), t1(ev), somaCtx.spatialSize, somaCtx.lmkMasks{lmk}, lmkCenter);
            mp.soma_centroid_dist_px = somaMetric.soma_centroid_dist_px;
            mp.soma_start_dist_px = somaMetric.soma_start_dist_px;
            mp.soma_end_dist_px = somaMetric.soma_end_dist_px;
            mp.landmark_dist_avg = matrixAt(somaCtx.distAvg, ev, lmk);
            mp.landmark_dist_min = matrixAt(somaCtx.distMin, ev, lmk);
            mp.chg_toward = matrixAt(somaCtx.chgToward, ev, lmk);
            mp.chg_away = matrixAt(somaCtx.chgAway, ev, lmk);
            mp.chg_toward_bef_reach = matrixAt(somaCtx.chgTowardBefReach, ev, lmk);
            mp.chg_away_aft_reach = matrixAt(somaCtx.chgAwayAftReach, ev, lmk);
            den = mp.chg_toward + mp.chg_away;
            if isfinite(den) && den > 0
                mp.motion_index = (mp.chg_toward - mp.chg_away) / den;
            end
            mp.motion_class = classifyMotion(mp.chg_toward, mp.chg_away, mp.chg_toward_bef_reach, mp.chg_away_aft_reach);
        else
            mp.motion_class = "no_landmark";
        end
        eventMapRows(rr,1) = mp;
    end
    for ee = 1:nEvents
        idxEvt = find(mapPairs(:,1) == ee);
        if isempty(idxEvt)
            continue
        end
        primaryHit = idxEvt([eventMapRows(idxEvt).is_primary_cell]);
        if isempty(primaryHit)
            continue
        end
        primaryEntry = eventMapRows(primaryHit(1)).cell_first_overlap_frame;
        primaryEntryThresh = eventMapRows(primaryHit(1)).cell_entry_frame_cov_thresh;
        if ~isfinite(primaryEntry)
            continue
        end
        for jj = idxEvt(:)'
            if isfinite(eventMapRows(jj).cell_first_overlap_frame)
                eventMapRows(jj).cell_entry_delay_from_primary_frames = ...
                    eventMapRows(jj).cell_first_overlap_frame - primaryEntry;
            end
        end
        if ~isfinite(primaryEntryThresh)
            primaryEntryThresh = primaryEntry;
        end
        if isfinite(primaryEntryThresh)
            for jj = idxEvt(:)'
                cellEntryThresh = eventMapRows(jj).cell_entry_frame_cov_thresh;
                if isfinite(cellEntryThresh)
                    eventMapRows(jj).cell_entry_delay_from_primary_cov_thresh_frames = ...
                        cellEntryThresh - primaryEntryThresh;
                elseif eventMapRows(jj).is_primary_cell
                    eventMapRows(jj).cell_entry_delay_from_primary_cov_thresh_frames = 0;
                end
            end
        end
    end
    mRows = localFilterOfficialMapRows(eventMapRows, minAssociatedCellCoverageFrac);
end

function cls = classifyMotion(toward, away, towardBeforeReach, awayAfterReach)
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

function v = matrixAt(x, i, j)
    v = NaN;
    if isempty(x) || ~isnumeric(x) || i < 1 || j < 1
        return
    end
    if size(x,1) >= i && size(x,2) >= j
        v = x(i,j);
    end
end

function p = eventSpatialPixels(ctx, ev)
    p = [];
    if isfield(ctx, 'eventPix') && iscell(ctx.eventPix) && numel(ctx.eventPix) >= ev
        p = ctx.eventPix{ev};
        if ~isempty(p)
            return
        end
    end
    if ev < 1 || ev > numel(ctx.xSpa) || ev > numel(ctx.xSpaTemp)
        return
    end
    if iscell(ctx.xSpa) && numel(ctx.xSpa) >= ev
        px = ctx.xSpa{ev};
        if isnumeric(px)
            p = unique(px(:));
            p = p(p > 0 & isfinite(p));
            if ~isempty(p)
                return
            end
        end
    end
    if iscell(ctx.xSpaTemp) && numel(ctx.xSpaTemp) >= ev
        px = ctx.xSpaTemp{ev};
        if isnumeric(px)
            p = voxelToSpatial(px, ctx.spatialSize);
        end
    end
end

function p = voxelToSpatial(voxelIdx, spatialSize)
    p = [];
    if isempty(voxelIdx) || numel(spatialSize) < 3 || any(~isfinite(spatialSize(1:3)))
        return
    end
    sx = round(spatialSize(1));
    sy = round(spatialSize(2));
    sz = round(spatialSize(3));
    if sx <= 0 || sy <= 0 || sz <= 0
        return
    end
    try
        [ih, iw, il, ~] = ind2sub([sx, sy, sz, max(1,round(spatialSize(4)))], voxelIdx(:));
        p = unique(sub2ind([sx, sy, sz], ih, iw, il));
    catch
        p = [];
    end
end

function ctx = buildSomaContext(fts1, nEvents, spatialSizeHint)
    ctx = struct();
    ctx.cellToSoma = NaN(0,1);
    ctx.lmkPix = {};
    ctx.lmkMasks = {};
    ctx.cellMasks = {};
    ctx.xSpa = getPathOrDefault(fts1, 'loc.xSpa', {});
    ctx.xSpaTemp = getPathOrDefault(fts1, 'loc.xSpaTemp', {});
    ctx.spatialSize = [NaN NaN NaN NaN];
    if nargin >= 3 && ~isempty(spatialSizeHint)
        spatialSizeHint = double(spatialSizeHint);
        ctx.spatialSize(1:min(numel(spatialSizeHint), 4)) = spatialSizeHint(1:min(numel(spatialSizeHint), 4));
    end
    ctx.distAvg = force2d(getPathOrDefault(fts1, 'region.landmarkDist.distAvg', []), nEvents);
    ctx.distMin = force2d(getPathOrDefault(fts1, 'region.landmarkDist.distMin', []), nEvents);
    ctx.chgToward = force2d(getPathOrDefault(fts1, 'region.landmarkDir.chgToward', []), nEvents);
    ctx.chgAway = force2d(getPathOrDefault(fts1, 'region.landmarkDir.chgAway', []), nEvents);
    ctx.chgTowardBefReach = force2d(getPathOrDefault(fts1, 'region.landmarkDir.chgTowardBefReach', []), nEvents);
    ctx.chgAwayAftReach = force2d(getPathOrDefault(fts1, 'region.landmarkDir.chgAwayAftReach', []), nEvents);

    masks = getPathOrDefault(fts1, 'region.landMark.mask', {});
    cellMasks = getPathOrDefault(fts1, 'region.cell.mask', {});
    ctx.cellMasks = cellMasks;
    ctx.cellPix = {};
    inclu = getPathOrDefault(fts1, 'region.cell.incluLmk', []);
    cellCenters = getPathOrDefault(fts1, 'region.cell.center', []);
    lmkCenters = getPathOrDefault(fts1, 'region.landMark.center', []);
    ctx.lmkCenters = lmkCenters;
    if iscell(masks)
        ctx.lmkMasks = masks;
        ctx.lmkPix = cell(numel(masks), 1);
        for ii = 1:numel(masks)
            m = masks{ii};
            if isnumeric(m) || islogical(m)
                ctx.lmkPix{ii} = find(m > 0);
                if numel(ctx.spatialSize) >= 3 && any(~isfinite(ctx.spatialSize(1:3)))
                    sz = size(m);
                    ctx.spatialSize(1:min(3,numel(sz))) = sz(1:min(3,numel(sz)));
                end
            else
                ctx.lmkPix{ii} = [];
            end
        end
    end
    if iscell(cellMasks)
        ctx.cellPix = cell(numel(cellMasks), 1);
        for ii = 1:numel(cellMasks)
            m = cellMasks{ii};
            if isnumeric(m) || islogical(m)
                ctx.cellPix{ii} = find(m > 0);
                if numel(ctx.spatialSize) >= 3 && any(~isfinite(ctx.spatialSize(1:3)))
                    sz = size(m);
                    ctx.spatialSize(1:min(3,numel(sz))) = sz(1:min(3,numel(sz)));
                end
            else
                ctx.cellPix{ii} = [];
            end
        end
    end

    ctx.eventPix = cell(nEvents, 1);
    for ee = 1:nEvents
        if ee < 1 || ee > numel(ctx.xSpa) || ee > numel(ctx.xSpaTemp)
            ctx.eventPix{ee} = [];
            continue
        end
        if iscell(ctx.xSpa) && numel(ctx.xSpa) >= ee
            px = ctx.xSpa{ee};
            if isnumeric(px)
                p = unique(px(:));
                p = p(p > 0 & isfinite(p));
                if ~isempty(p)
                    ctx.eventPix{ee} = p;
                    continue
                end
            end
        end
        if iscell(ctx.xSpaTemp) && numel(ctx.xSpaTemp) >= ee
            px = ctx.xSpaTemp{ee};
            if isnumeric(px)
                ctx.eventPix{ee} = voxelToSpatial(px, ctx.spatialSize);
            else
                ctx.eventPix{ee} = [];
            end
        else
            ctx.eventPix{ee} = [];
        end
    end

    ctx.cellToSoma = resolveSingleSomaPerCell(inclu, cellMasks, masks, cellCenters, lmkCenters);
end

function out = safeCellFetch(C, idx)
out = [];
if iscell(C) && numel(C) >= idx
    out = C{idx};
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
        bestLmk = pickMostCentralLandmark(cand, cc, cellMasks, lmkMasks, cellCenters, lmkCenters);
        cellToSoma(cc) = bestLmk;
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
                scoreOverlap(ii) = overlapCountSafe(cm, lm);
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

function n = overlapCountSafe(a, b)
    n = -inf;
    if ~(isnumeric(a) || islogical(a)) || ~(isnumeric(b) || islogical(b))
        return
    end
    if isequal(size(a), size(b))
        n = nnz((a > 0) & (b > 0));
        return
    end
    ai = find(a > 0);
    bi = find(b > 0);
    if isempty(ai) || isempty(bi)
        n = 0;
        return
    end
    n = numel(intersect(ai(:), bi(:)));
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

function out = force2d(x, nRows)
    if isempty(x) || ~isnumeric(x)
        out = NaN(nRows, 0);
        return
    end
    out = x;
    if isvector(out)
        out = out(:);
    end
    if size(out,1) < nRows
        tmp = NaN(nRows, size(out,2));
        tmp(1:size(out,1),:) = out;
        out = tmp;
    end
end

function [pairs, primaryCellID] = eventCellPairs(fts1, nEvents)
    primaryCellID = NaN(nEvents,1);
    pairs = zeros(0,3);
    memberIdx = getPathOrDefault(fts1, 'region.cell.memberIdx', []);
    if isempty(memberIdx)
        return
    end

    if isnumeric(memberIdx)
        if isvector(memberIdx) && numel(memberIdx) == nEvents
            ids = forceColumn(memberIdx);
            primaryCellID = ids;
            valid = ~isnan(ids) & ids > 0;
            evt = find(valid);
            cids = ids(valid);
            pairs = [evt, cids, nan(numel(evt),1)];
            return
        end
        if size(memberIdx,1) == nEvents
            for ee = 1:nEvents
                row = memberIdx(ee,:);
                valid = find(row > 0);
                if isempty(valid)
                    continue
                end
                primaryCellID(ee) = valid(1);
                scores = row(valid)';
                pairs = [pairs; [repmat(ee, numel(valid), 1), valid', scores]]; %#ok<AGROW>
            end
            return
        end
    end

    if iscell(memberIdx) && numel(memberIdx) == nEvents
        for ee = 1:nEvents
            ids = memberIdx{ee};
            if isempty(ids)
                continue
            end
            ids = ids(:);
            ids = ids(~isnan(ids) & ids > 0);
            if isempty(ids)
                continue
            end
            primaryCellID(ee) = ids(1);
            pairs = [pairs; [repmat(ee, numel(ids), 1), ids, nan(numel(ids),1)]]; %#ok<AGROW>
        end
    end
end

function pairs = localEnsurePrimaryPairs(pairs, primaryCellID)
if isempty(primaryCellID)
    return
end
for ee = 1:numel(primaryCellID)
    cc = primaryCellID(ee);
    if ~(isfinite(cc) && cc > 0)
        continue
    end
    if isempty(pairs)
        pairs = [ee, cc, NaN]; %#ok<AGROW>
        continue
    end
    mask = pairs(:,1) == ee & pairs(:,2) == cc;
    if any(mask)
        continue
    end
    pairs = [pairs; ee, cc, NaN]; %#ok<AGROW>
end
end

function rows = localFilterOfficialMapRows(rows, minCoverageFrac)
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
rows = rows(isPrimary | (isfinite(peakFrac) & peakFrac >= minCoverageFrac));
end

function ownerCell = localResolveOwnerCellFromPointSingle(pt, cellMasks, candidateCells)
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
    ownerCell = localNearestCellMaskCentroidSingle(pt, candidateCells, cellMasks);
else
    ownerCell = localNearestCellMaskCentroidSingle(pt, inside, cellMasks);
end
end

function ownerCell = localNearestCellMaskCentroidSingle(pt, candidateCells, cellMasks)
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

function out = getPathOrDefault(s, dottedPath, defaultVal)
    out = defaultVal;
    if ~isstruct(s) || isempty(dottedPath)
        return
    end
    parts = strsplit(dottedPath, '.');
    cur = s;
    for ii = 1:numel(parts)
        p = parts{ii};
        if ~isstruct(cur) || ~isfield(cur, p)
            return
        end
        cur = cur.(p);
    end
    out = cur;
end

function v = forceColumn(x)
    if isempty(x)
        v = [];
    else
        v = x(:);
    end
end

function v = forceColumnWithLength(x, n)
    if isempty(x)
        v = NaN(n,1);
        return
    end
    xv = x(:);
    if numel(xv) == n
        v = xv;
    else
        v = NaN(n,1);
        m = min(numel(xv), n);
        v(1:m) = xv(1:m);
    end
end

function out = scalarOrNaN(x)
    if isempty(x)
        out = NaN;
        return
    end
    if isnumeric(x)
        out = x(1);
    else
        out = NaN;
    end
end

function s = asString(v)
    if isstring(v)
        s = v;
    elseif ischar(v)
        s = string(v);
    elseif isempty(v)
        s = "";
    else
        s = string(v);
    end
end

function c = asChar(v)
    if ischar(v)
        c = v;
    elseif isstring(v)
        c = char(v);
    elseif isempty(v)
        c = '';
    else
        c = char(string(v));
    end
end

function out = tf(v)
    if v
        out = 'true';
    else
        out = 'false';
    end
end
