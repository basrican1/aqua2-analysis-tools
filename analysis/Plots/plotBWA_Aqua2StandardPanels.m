function out = plotBWA_Aqua2StandardPanels(ctx, varargin)
% Render standard legacy BWA_Aqua2 panel set from precomputed context.
% Input:
%   ctx: output struct from prepareLegacyAqua2Context()

p = inputParser;
p.StructExpand = false;
addRequired(p, 'ctx', @isstruct);
addParameter(p, 'style', struct(), @isstruct);
addParameter(p, 'statsConfig', struct(), @isstruct);
parse(p, ctx, varargin{:});
style = p.Results.style;
statsConfig = p.Results.statsConfig;
ctx = p.Results.ctx;

    if ~isstruct(ctx) || ~isfield(ctx, 'BWA_AQ2') || ~isfield(ctx, 'opts')
        error('plotBWA_Aqua2StandardPanels:BadContext', ...
            'ctx must be output from prepareLegacyAqua2Context().');
    end

    BWA_AQ2 = ctx.BWA_AQ2;
    opts = ctx.opts;
    timeFrameRangeEpochs = ctx.timeFrameRangeEpochs;
    memberlist = ctx.memberlistRow3;
    associatedCoverageTraces = ctx.associatedCoverageTraces;
    memberCoverage = ctx.memberCoverage;
    filesTbl = table();
    eventsTbl = table();
    if isfield(ctx, 'filesTbl'), filesTbl = ctx.filesTbl; end
    if isfield(ctx, 'eventsTbl'), eventsTbl = ctx.eventsTbl; end

    out = struct();
    out.basic.area = BAviolins(BWA_AQ2.fts1, opts, 'basic.area', timeFrameRangeEpochs, memberlist, ...
        'style', style, 'statsConfig', statsConfig, 'eventsTbl', eventsTbl, 'filesTbl', filesTbl);
    out.curve.duration = BAviolins(BWA_AQ2.fts1, opts, 'curve.duration', timeFrameRangeEpochs, memberlist, ...
        'style', style, 'statsConfig', statsConfig, 'eventsTbl', eventsTbl, 'filesTbl', filesTbl);
    out.curve.rise19 = BAviolins(BWA_AQ2.fts1, opts, 'curve.rise19', timeFrameRangeEpochs, memberlist, ...
        'style', style, 'statsConfig', statsConfig, 'eventsTbl', eventsTbl, 'filesTbl', filesTbl);
    out.curve.fall91 = BAviolins(BWA_AQ2.fts1, opts, 'curve.fall91', timeFrameRangeEpochs, memberlist, ...
        'style', style, 'statsConfig', statsConfig, 'eventsTbl', eventsTbl, 'filesTbl', filesTbl);
    out.curve.decayTau = BAviolins(BWA_AQ2.fts1, opts, 'curve.decayTau', timeFrameRangeEpochs, memberlist, ...
        'style', style, 'statsConfig', statsConfig, 'eventsTbl', eventsTbl, 'filesTbl', filesTbl);
    out.curve.dffAUC = BAviolins(BWA_AQ2.fts1, opts, 'curve.dffAUC', timeFrameRangeEpochs, memberlist, ...
        'style', style, 'statsConfig', statsConfig, 'eventsTbl', eventsTbl, 'filesTbl', filesTbl);
    out.curve.dffMax2 = BAviolins(BWA_AQ2.fts1, opts, 'curve.dffMax2', timeFrameRangeEpochs, memberlist, ...
        'style', style, 'statsConfig', statsConfig, 'eventsTbl', eventsTbl, 'filesTbl', filesTbl);

    doBySlice = 0;
    areaPercentiles = [0 100];
    out.freq = nEventsPerMin(BWA_AQ2.fts1, opts, memberlist, timeFrameRangeEpochs, areaPercentiles, doBySlice);
    covThreshold = 0.1;
    [coverage, t0_subset, epMembIdx] = CellCoverage_calc(BWA_AQ2.fts1, opts, timeFrameRangeEpochs, memberlist, covThreshold, associatedCoverageTraces);
    tOccuring_subset = t0_subset * opts.frameRate(1) / 60;
    out.coverage = CellCoverage_plot(coverage, tOccuring_subset, covThreshold, timeFrameRangeEpochs, epMembIdx);
    out.footprint = pCellFoot(BWA_AQ2.fts1, opts, timeFrameRangeEpochs, memberlist, memberCoverage);

    [out.propagationDistance, out.propagationAngle, out.propagationDistanceFigure, out.propagationAngleFigure] = ...
        propagation(BWA_AQ2.fts1, opts, timeFrameRangeEpochs, memberlist);
    out.propagationSpeed = BAviolins(BWA_AQ2.fts1, opts, 'propagation.maxPropSpeed', timeFrameRangeEpochs, memberlist, ...
        'style', style, 'statsConfig', statsConfig, 'eventsTbl', eventsTbl, 'filesTbl', filesTbl);

    allCentroids = CalcCentroids(BWA_AQ2.fts1.loc.xSpaTemp, BWA_AQ2.fts1.loc.t0, BWA_AQ2.fts1.loc.t1, opts.sz(1,:), BWA_AQ2.fts1.propagation.propGrow);
    [out.centroidDist, out.centroidJumps, out.centroidDistTot] = DistFromCentroids(allCentroids, opts.spatialRes(1));
    out.centroidFigure = distCentroidViolin(BWA_AQ2.fts1, opts, out.centroidDistTot, timeFrameRangeEpochs, memberlist);
end
