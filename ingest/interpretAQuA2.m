function [BWA_AQ2, canonicalOut] = interpretAQuA2(varargin)
    % Main entrypoint: always run legacy BWA_AQ2 accumulation + canonical ingest.
    %
    % Usage:
    %   1) Default (legacy + canonical):
    %      [BWA_AQ2, canonicalOut] = interpretAQuA2();
    %
    %   2) With parallel legacy extraction (cluster multi-CPU):
    %      [BWA_AQ2, canonicalOut] = interpretAQuA2('useParallelLegacy', true);
    %
    %   3) With parallel canonical ingest progress:
    %      [BWA_AQ2, canonicalOut] = interpretAQuA2( ...
    %          'useParallelLegacy', true, 'useParallelCanonical', true);
    %
    %   4) Custom canonical output location and CSV export:
    %      [BWA_AQ2, canonicalOut] = interpretAQuA2( ...
    %          'canonicalOutputRoot', fullfile('derived', 'ingest_runs'), ...
    %          'canonicalSaveCsv', false);
    %
    % Notes:
    %   - canonicalOut.runDir points to the timestamped ingest folder.
    %   - 'canonicalSaveCsv' writes CSVs in addition to MAT tables.
    %   - Canonical is always built from in-memory BWA_AQ2 to avoid rereading AQuA2 files.
    %   - For local analysis, copy the canonical run folder to local:
    %       files.mat, cells.mat, events.mat, event_cell_map.mat, manifest.mat/json
    %   - event_cell_map now includes soma/landmark metrics when available:
    %       soma overlap fractions, soma>=50% flag, toward/away motion fields,
    %       motion_index and motion_class, plus cell-overlap / outside-cell fractions.
    %   - event_cell_map also carries event-cell partition timing metrics:
    %       first/last overlap frame, overlap duration, peak overlap, and
    %       entry delay from the primary/origin cell.
    %   - events now include perimeter/circularity and continuous propagation-distance metrics.
    %   - primary_cell_uid / is_primary_cell now mean onset-based origin cell among
    %     the AQuA2 region.cell.memberIdx candidate cells, not "first positive memberIdx".
    %     Optional: *.csv (if canonicalSaveCsv=true), latest_ingest.txt (pointer only).
    %   - Pipeline calls CreateTopDataBases with saveLegacyDbs=false, so
    %     DBfiles.mat / DBanimals.mat / DBExpCond.mat are not rewritten.
    %   - BWA_AQ2.mat / DBfiles.mat / DBanimals.mat are not required for canonical-driven local analysis.
    BWA_setupPaths(fileparts(mfilename('fullpath')));

    p = inputParser;
    p.StructExpand = false;
    addParameter(p, 'canonicalOutputRoot', fullfile('derived', 'ingest_runs'), @ischar);
    addParameter(p, 'canonicalSaveCsv', false, @islogical);
    addParameter(p, 'useParallelLegacy', false, @islogical);
    addParameter(p, 'useParallelCanonical', false, @islogical);
    parse(p, varargin{:});

    canonicalOutputRoot = p.Results.canonicalOutputRoot;
    canonicalSaveCsv = p.Results.canonicalSaveCsv;
    useParallelLegacy = p.Results.useParallelLegacy;
    useParallelCanonical = p.Results.useParallelCanonical;

    disp('interpretAQuA2: starting run');
    disp(['  useParallelLegacy=' tf(useParallelLegacy) ...
        ', useParallelCanonical=' tf(useParallelCanonical) ...
        ', canonicalOutputRoot=' canonicalOutputRoot ...
        ', canonicalSaveCsv=' tf(canonicalSaveCsv)]);

    BWA_AQ2 = [];
    canonicalOut = [];
    DB_Files = [];
    DB_Animals = [];
    DB_ExpCond = struct('fileID',{},'animalID',{},'regionImaged',{},...
        'virus',{},'virusLoc',{},'dpi_virus',{},...
        'drug',{}, 'drugLoc', {}, 'dpi_drug',{},...
        'stim', {}, 'stimloc', {}, 'dpi_stim', {},...
        'training', {}, 'dpi_training', {},...
        'imageEpoc', {}, 'epocFrameStart', {}, 'epocFrameEnd', {}, 'epocTimeStart', {}, 'epocTimeEnd',{});

    %recursive search, make database lists
    tLegacy = tic;
    disp('Legacy: building databases...');
    [DB_Files, DB_Animals, DB_ExpCond] = CreateTopDataBases('saveLegacyDbs', false);

    %Create empty Structure
    BWA_AQ2 = MakeAQ2();
    fieldList = getAllFieldsAtLevel(BWA_AQ2, 5);

    %loop and summarize data from all entries
    nItems = size(DB_Files,2);
    disp(['Legacy: processing ' num2str(nItems) ' files']);
    newDatCell = cell(nItems,1);
    if useParallelLegacy
        didParallel = false;
        try
            if isempty(gcp('nocreate'))
                parpool;
            end
            pObj = gcp('nocreate');
            disp(['Legacy: using parallel pool with ' num2str(pObj.NumWorkers) ' workers']);
            q = parallel.pool.DataQueue;
            nDone = 0;
            reportStep = max(1, floor(nItems/20));
            afterEach(q, @legacyProgressTick);
            parfor ii=1:nItems
                thisFileInfo = DB_Files(ii);
                IDStruct = getIDfileBy(thisFileInfo.filepath, 'fileID', thisFileInfo.fileID);
                newDatCell{ii} = extractFromAQ2(thisFileInfo, IDStruct, fieldList);
                send(q, ii);
            end
            didParallel = true;
        catch
            didParallel = false;
        end
        if ~didParallel
            disp('Legacy: parallel unavailable, falling back to sequential mode');
            for ii=1:nItems
                thisFileInfo = DB_Files(ii);
                disp(['Legacy seq [' num2str(ii) '/' num2str(nItems) ']: ' thisFileInfo.filename]);
                IDStruct = getIDfileBy(thisFileInfo.filepath, 'fileID', thisFileInfo.fileID);
                newDatCell{ii} = extractFromAQ2(thisFileInfo, IDStruct, fieldList);
            end
        end
        disp('Legacy: merging extracted file structures...');
        for ii=1:nItems
            BWA_AQ2 = insertDatToSummary(BWA_AQ2, newDatCell{ii}, 5);
        end
    else
        for ii=1:nItems
            thisFileInfo = DB_Files(ii);
            disp(['Legacy seq [' num2str(ii) '/' num2str(nItems) ']: ' thisFileInfo.filename]);
            IDStruct = getIDfileBy(thisFileInfo.filepath, 'fileID', thisFileInfo.fileID);
            newDat = extractFromAQ2(thisFileInfo, IDStruct, fieldList);
            BWA_AQ2 = insertDatToSummary(BWA_AQ2, newDat, 5);
        end
    end
    disp(['Legacy: done in ' num2str(toc(tLegacy), '%.1f') ' s']);

    tCanonical = tic;
    disp('Canonical: building from in-memory BWA_AQ2 (no raw reread)');
    canonicalOut = ingestCanonicalFromBWA(BWA_AQ2, DB_Files, DB_Animals, DB_ExpCond, ...
        'outputRoot', canonicalOutputRoot, 'saveCsv', canonicalSaveCsv, ...
        'useParallel', useParallelCanonical);
    disp(['Canonical: done in ' num2str(toc(tCanonical), '%.1f') ' s']);

    % Save a slim legacy compatibility artifact after canonical ownership is resolved.
    BWA_AQ2_disk = BWA_AQ2;
    BWA_AQ2_disk = augmentLegacyOwnershipFromCanonical(BWA_AQ2_disk, canonicalOut);
    if isfield(BWA_AQ2_disk, 'DB')
        BWA_AQ2_disk = rmfield(BWA_AQ2_disk, 'DB');
    end
    if isfield(BWA_AQ2_disk, 'opts')
        BWA_AQ2_disk = rmfield(BWA_AQ2_disk, 'opts');
    end
    if isfield(BWA_AQ2_disk, 'fts1')
        if isfield(BWA_AQ2_disk.fts1, 'fileID')
            BWA_AQ2_disk.fts1 = rmfield(BWA_AQ2_disk.fts1, 'fileID');
        end
        if isfield(BWA_AQ2_disk.fts1, 'loc')
            if isfield(BWA_AQ2_disk.fts1.loc, 't0')
                BWA_AQ2_disk.fts1.loc = rmfield(BWA_AQ2_disk.fts1.loc, 't0');
            end
            if isfield(BWA_AQ2_disk.fts1.loc, 't1')
                BWA_AQ2_disk.fts1.loc = rmfield(BWA_AQ2_disk.fts1.loc, 't1');
            end
        end
    end
    % Strip event metric fields already represented in canonical events table.
    % Keep structural fields (xSpaTemp/memberIdx/masks/propGrow/propShrink) for legacy-style analyses.
    BWA_AQ2_disk = stripCanonicalBackedMetrics(BWA_AQ2_disk);
    tmpSave = struct();
    tmpSave.BWA_AQ2 = BWA_AQ2_disk;
    save('BWA_AQ2.mat', '-struct', 'tmpSave');

    disp('interpretAQuA2: finished');

    function legacyProgressTick(~)
        nDone = nDone + 1;
        if nDone == 1 || nDone == nItems || mod(nDone, reportStep) == 0
            disp(['Legacy par [' num2str(nDone) '/' num2str(nItems) ']']);
        end
    end
end

function out = tf(v)
    if v
        out = 'true';
    else
        out = 'false';
    end
end



function structOut = MakeAQ2()
    % Lean legacy payload: keep only fields not fully captured in canonical tables.
    % NOTE:
    %   Canonical ingest is built from in-memory BWA_AQ2 via ingestCanonicalFromBWA().
    %   Therefore, event metric fields must still be present here so extractFromAQ2()
    %   captures them from raw res. If omitted, canonical event metrics become NaN.
    structOut.DB.fileID = '';
    structOut.DB.animalID = '';
    structOut.DB.sliceNum = '';
    structOut.DB.nEvents = [];
    structOut.DB.nCells = [];
    structOut.DB.nSomas = [];
    
    structOut.fts1.loc.t0 = [];
    structOut.fts1.loc.t1 = [];
    structOut.fts1.loc.xSpa = [];
    structOut.fts1.loc.xSpaTemp = [];
    structOut.fts1.basic.area = [];
    structOut.fts1.curve.duration = [];
    structOut.fts1.curve.rise19 = [];
    structOut.fts1.curve.fall91 = [];
    structOut.fts1.curve.decayTau = [];
    structOut.fts1.curve.dffAUC = [];
    structOut.fts1.curve.dffMax2 = [];
    structOut.fts1.region.cell.memberIdx = [];
    structOut.fts1.region.cell.mask = [];
    structOut.fts1.region.cell.incluLmk = [];
    structOut.fts1.region.landMark.mask = [];
    structOut.fts1.region.landmarkDist.distAvg = [];
    structOut.fts1.region.landmarkDist.distMin = [];
    structOut.fts1.region.landmarkDir.chgToward = [];
    structOut.fts1.region.landmarkDir.chgAway = [];
    structOut.fts1.region.landmarkDir.chgTowardBefReach = [];
    structOut.fts1.region.landmarkDir.chgAwayAftReach = [];
    structOut.fts1.propagation.maxPropSpeed = [];
    structOut.fts1.propagation.propGrow = [];
    structOut.fts1.propagation.propShrink = [];

    structOut.opts.sz = [];
    structOut.opts.frameRate = [];
    structOut.opts.spatialRes = [];

end

function S = stripCanonicalBackedMetrics(S)
    % Remove metric branches duplicated in canonical events table.
    % Keep structural branches needed for legacy-style structural workflows.
    paths = { ...
        'fts1.basic', ...
        'fts1.curve', ...
        'fts1.propagation.maxPropSpeed'};
    for ii = 1:numel(paths)
        S = rmFieldPathIfExists(S, paths{ii});
    end
end

function B = augmentLegacyOwnershipFromCanonical(B, canonicalOut)
if ~(isstruct(B) && isfield(B, 'fts1') && isstruct(B.fts1) && isstruct(canonicalOut))
    return
end
E = [];
if isfield(canonicalOut, 'events') && istable(canonicalOut.events)
    E = canonicalOut.events;
elseif isfield(canonicalOut, 'subset') && isstruct(canonicalOut.subset) && ...
        isfield(canonicalOut.subset, 'events') && istable(canonicalOut.subset.events)
    E = canonicalOut.subset.events;
end
if isempty(E)
    return
end
need = {'event_uid','primary_cell_uid','branch_origin_x_px','branch_origin_y_px','branch_origin_frame','branch_origin_cell_uid'};
if ~all(ismember(need, E.Properties.VariableNames))
    return
end
if ~isfield(B.fts1, 'fileID') || ~isfield(B.fts1, 'loc') || ~isfield(B.fts1.loc, 'xSpaTemp')
    return
end
nEvents = numel(B.fts1.loc.xSpaTemp);
if nEvents < 1
    return
end
eventUID = strings(nEvents,1);
fileIDs = localCanonCellstrInterpret(B.fts1.fileID);
counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ii = 1:nEvents
    fid = fileIDs{ii};
    if ~isKey(counts, fid)
        counts(fid) = 0;
    end
    counts(fid) = counts(fid) + 1;
    eventUID(ii) = string(sprintf('%s_E%d', fid, counts(fid)));
end
lookup = containers.Map('KeyType', 'char', 'ValueType', 'double');
canonUID = string(E.event_uid);
for ii = 1:height(E)
    lookup(char(canonUID(ii))) = ii;
end
primaryUID = repmat("", nEvents, 1);
branchOwnerUID = repmat("", nEvents, 1);
branchFrame = NaN(nEvents, 1);
branchXY = NaN(nEvents, 2);
for ii = 1:nEvents
    key = char(eventUID(ii));
    if ~isKey(lookup, key)
        continue
    end
    jj = lookup(key);
    primaryUID(ii) = string(E.primary_cell_uid(jj));
    branchOwnerUID(ii) = string(E.branch_origin_cell_uid(jj));
    branchFrame(ii) = double(E.branch_origin_frame(jj));
    branchXY(ii,1) = double(E.branch_origin_x_px(jj));
    branchXY(ii,2) = double(E.branch_origin_y_px(jj));
end
B.fts1.primary_cell_uid = primaryUID;
B.fts1.branch_origin_cell_uid = branchOwnerUID;
if ~isfield(B.fts1, 'branch_origin')
    B.fts1.branch_origin = struct();
end
B.fts1.branch_origin.frame = branchFrame;
B.fts1.branch_origin.xy = branchXY;
end

function c = localCanonCellstrInterpret(v)
if iscell(v)
    c = cell(size(v));
    for i = 1:numel(v)
        c{i} = char(string(v{i}));
    end
elseif isstring(v)
    c = cellstr(v(:));
elseif ischar(v)
    c = {v};
else
    c = cellstr(string(v(:)));
end
end

function S = rmFieldPathIfExists(S, dottedPath)
    parts = strsplit(dottedPath, '.');
    if isempty(parts)
        return
    end
    S = rmFieldParts(S, parts);
end

function out = rmFieldParts(in, parts)
    out = in;
    if ~isstruct(out) || isempty(parts)
        return
    end
    f = parts{1};
    if ~isfield(out, f)
        return
    end
    if numel(parts) == 1
        out = rmfield(out, f);
        return
    end
    child = out.(f);
    child = rmFieldParts(child, parts(2:end));
    out.(f) = child;
end
