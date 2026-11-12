function [h_fig, means_alln, means_celln] = plotBAviolinPanels(comp, timeFrameRangeEpochs, memberlist, varargin)
% Plot-only stage for BAviolins using precomputed comp struct.
% This function is the large organizer. Individual subplot functions can
% also be called standalone.
%
% The top-right inset slice spider can either use a fixed epoch pair
% (legacy behavior) or a combined post-epoch definition passed through
% spiderCombineEpochNames/spiderCombinedEpochLabel.

p = inputParser;
p.StructExpand = false;
addRequired(p, 'comp', @isstruct);
addRequired(p, 'timeFrameRangeEpochs', @isstruct);
addRequired(p, 'memberlist', @isnumeric);
addParameter(p,'optFigSuffix', '', @ischar);
addParameter(p,'doByCellOnly', false, @islogical);
addParameter(p,'style', struct(), @isstruct);
addParameter(p,'statsConfig', struct(), @isstruct);
addParameter(p,'showSliceSpider', true, @islogical);
addParameter(p,'spiderBranchMode', 'both', @(x) any(strcmpi(string(x), ["cell_only","both"])));
addParameter(p,'spiderCombineEpochNames', {}, @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p,'spiderCombinedEpochLabel', 'CNO_post', @(x) ischar(x) || isstring(x));
parse(p, comp, timeFrameRangeEpochs, memberlist, varargin{:});
style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);
statsConfig = p.Results.statsConfig;
doByCellOnly = p.Results.doByCellOnly;
showSliceSpider = p.Results.showSliceSpider;
spiderBranchMode = char(string(p.Results.spiderBranchMode));
if doByCellOnly
    means_alln = NaN;
else
    means_alln = [];
end

if isempty(p.Results.optFigSuffix)
    optFigSuffix = '';
else
    optFigSuffix = strcat('_', bwaFormatFileIDLabel(p.Results.optFigSuffix));
end

datTxt = comp.datTxt;

fName = strcat('BAVio_', datTxt, optFigSuffix);
fUniqueTag = 'BA_violins';
h_fig = findall(0, 'Type', 'figure', 'Tag', fUniqueTag, 'Name', fName);
if isempty(h_fig)
    h_fig = figure('Position', style.figurePositionBA, 'tag', fUniqueTag, 'Name', fName);
else
    figure(h_fig);
end

ax = createPanelAxes(h_fig);

if ~doByCellOnly
    plotBAviolinSubplotTimeScatter(comp.dffAll, comp.tOccuring, memberlist, datTxt, timeFrameRangeEpochs, ...
        'figure', h_fig, 'axes', ax.ax1, 'style', style);
    ep = plotBAviolinSubplotEpochViolin(comp.dffAll, comp.t0, datTxt, timeFrameRangeEpochs, ...
        'figure', h_fig, 'axes', ax.ax2, 'style', style); %#ok<NASGU>

    tEpochs = timeFrameRangeEpochs.epochNames(comp.epMembIdx(~comp.nonEpoch));
    anovaAll = plotBAviolinSubplotAnova(comp.dffAll(~comp.nonEpoch), tEpochs, ...
        'figure', h_fig, 'axesPlot', ax.ax4, 'axesTable', ax.ax8, 'style', style, ...
        'statsConfig', statsConfig);
    mcAll = plotBAviolinSubplotMultcompare(anovaAll.stats, 'figure', h_fig, 'axes', ax.ax6, 'style', style);
    c_e = mcAll.c;
    m_e = mcAll.m;
else
    c_e = [];
    m_e = [];
end

figure(h_fig); axes(ax.ax3); cla(ax.ax3);
set(ax.ax3, 'FontSize', style.fontSizeAxis);
violinplot(comp.avgs_byCell, comp.tCat_byCell, 'ShowMedian', true, 'ShowMean', true, ...
    'GroupOrder', timeFrameRangeEpochs.epochNames);
ylabel(bwaMetricLabelLatex(datTxt, 'cell_raw'), 'Interpreter', 'latex');
anovaCell = plotBAviolinSubplotAnova(comp.avgs_byCell, comp.tCat_byCell, ...
    'figure', h_fig, 'axesPlot', ax.ax5, 'axesTable', ax.ax9, 'style', style, ...
    'statsConfig', statsConfig, 'pairID', comp.cells_byCell);
mcCell = plotBAviolinSubplotMultcompare(anovaCell.stats, 'figure', h_fig, 'axes', ax.ax7, 'style', style);
c_c = mcCell.c;
m_c = mcCell.m;

figure(h_fig);
linkaxes([ax.ax1, ax.ax2, ax.ax3], 'y');
hZ = zoom(h_fig); hP = pan(h_fig);
set(hZ, 'ActionPostCallback', @(src,evt) updateLabels(src,evt,ax.ax1));
set(hP, 'ActionPostCallback', @(src,evt) updateLabels(src,evt,ax.ax1));

if ~doByCellOnly
    plotBAviolinSubplotTimeScatter(comp.normDFFs, comp.tOccuringN, comp.outMember, ['norm ',datTxt], timeFrameRangeEpochs, ...
        'figure', h_fig, 'axes', ax.axn1, 'style', style);
    yline(0);
    ylabel(bwaMetricLabelLatex(datTxt, 'norm'), 'Interpreter', 'latex');

    epn = plotBAviolinSubplotEpochViolin(comp.normDFFs, comp.outT0s, ['norm ',datTxt], timeFrameRangeEpochs, ...
        'figure', h_fig, 'axes', ax.axn2, 'style', style); %#ok<NASGU>
    ylabel('${\sf log}_2{\sf FC}$', 'Interpreter', 'latex');

    maskNorm = ~isnan(comp.outEpochIndex);
    tEpochsn = timeFrameRangeEpochs.epochNames(comp.outEpochIndex(maskNorm));
    anovaAllN = plotBAviolinSubplotAnova(comp.normDFFs(maskNorm), tEpochsn, ...
        'figure', h_fig, 'axesPlot', ax.axn4, 'axesTable', ax.axn8, 'style', style, ...
        'statsConfig', statsConfig);
    means_alln = anovaAllN.stats.groupStatsM(:,2)';
    mcAllN = plotBAviolinSubplotMultcompare(anovaAllN.stats, 'figure', h_fig, 'axes', ax.axn6, 'style', style);
    c_en = mcAllN.c;
    m_en = mcAllN.m;
else
    c_en = [];
    m_en = [];
end

figure(h_fig); axes(ax.axn3); cla(ax.axn3);
set(ax.axn3, 'FontSize', style.fontSizeAxis);
violinplot(comp.avgsn_byCell, comp.tCatn_byCell, 'ShowMedian', true, 'ShowMean', true, ...
    'GroupOrder', timeFrameRangeEpochs.epochNames);
ylabel(bwaMetricLabelLatex(datTxt, 'cell_norm'), 'Interpreter', 'latex');
anovaCellN = plotBAviolinSubplotAnova(comp.avgsn_byCell, comp.tCatn_byCell, ...
    'figure', h_fig, 'axesPlot', ax.axn5, 'axesTable', ax.axn9, 'style', style, ...
    'statsConfig', statsConfig, 'pairID', comp.cellsn_byCell);
means_celln = anovaCellN.stats.groupStatsM(:,2)';
mcCellN = plotBAviolinSubplotMultcompare(anovaCellN.stats, 'figure', h_fig, 'axes', ax.axn7, 'style', style);
c_cn = mcCellN.c;
m_cn = mcCellN.m;

% Right column: top spider, middle summary, lower compact fold/value panels.
if showSliceSpider
    try
        if isfield(comp, 'normByEvent') && isfield(comp, 'epMembIdx') && isfield(comp, 'fileID')
            plotBAviolinSubplotSliceSpider(comp, timeFrameRangeEpochs, datTxt, ...
                'figure', h_fig, 'axes', ax.axsp, 'style', style, ...
                'branchMode', spiderBranchMode, ...
                'combineEpochNames', p.Results.spiderCombineEpochNames, ...
                'combinedEpochLabel', p.Results.spiderCombinedEpochLabel);
        else
            axes(ax.axsp); cla(ax.axsp); axis(ax.axsp, 'off');
        end
    catch ME
        warning('plotBAviolinPanels:SliceSpiderFailed', ...
            'Slice spider failed for %s: %s', datTxt, ME.message);
        axes(ax.axsp); cla(ax.axsp); axis(ax.axsp, 'off');
        text(ax.axsp, 0.02, 0.95, sprintf('Spider failed:\\n%s', ME.message), ...
            'Units', 'normalized', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top', 'Interpreter', 'none', 'FontSize', 8);
    end
else
    axes(ax.axsp); cla(ax.axsp); axis(ax.axsp, 'off');
end
plotBAviolinSubplotMetricSummary(comp, timeFrameRangeEpochs, datTxt, ...
    'figure', h_fig, 'axes', ax.axsum, 'style', style);

if ~doByCellOnly
    plotBAviolinSubplotFoldScatter(comp.outDFFs, comp.normDFFs, comp.outEpochIndex, datTxt, ...
        'figure', h_fig, 'axes', ax.axv1, 'style', style);
    plotBAviolinSubplotValueVsFold(c_e, m_e, c_en, m_en, datTxt, ...
        'figure', h_fig, 'axes', ax.axv3, 'style', style, ...
        'groupNames', mcAllN.groupNames);
else
    axes(ax.axv1); cla(ax.axv1); axis(ax.axv1, 'off');
    axes(ax.axv3); cla(ax.axv3); axis(ax.axv3, 'off');
end
plotBAviolinSubplotFoldScatter(comp.outDffs_c, comp.avgsn_byCell, comp.outEpochIndex_c, [datTxt '_c'], ...
    'figure', h_fig, 'axes', ax.axv2, 'style', style);
plotBAviolinSubplotValueVsFold(c_c, m_c, c_cn, m_cn, [datTxt '_c'], ...
    'figure', h_fig, 'axes', ax.axv4, 'style', style, ...
    'groupNames', mcCellN.groupNames);
end

function ax = createPanelAxes(h_fig)
buffer = 0.05;
bufferplus = 0.075;
nAcross = 6;
nbps = 2;
widthSingle = (1 - buffer - nAcross*buffer - nbps*(bufferplus-buffer))/nAcross;
widthDouble = 2*widthSingle+buffer;

nDown = 5;
heightSingle = (1 - buffer - nDown*buffer)/nDown;
plusminus = 0.02;
heightplus = heightSingle + plusminus;
heightminus = heightSingle - plusminus;

gL0 = buffer;
gL1 = gL0 +widthSingle + buffer;
gL2 = gL1 + widthSingle + bufferplus;
gL3 = gL2 + widthSingle + buffer;
gL4 = gL3 + widthSingle + bufferplus;
gL5 = gL4 + widthSingle + buffer;

gb0 = buffer;
gb1 = gb0 + heightSingle + buffer;
gb2 = gb1 + heightSingle + buffer;
gb3 = gb2 + heightminus + buffer;
gb4 = gb3 + heightSingle + buffer;

ax.ax1 = subplot('Position', [gL0, gb4, widthDouble, heightplus], 'Parent', h_fig);
ax.ax2 = subplot('Position', [gL0, gb3, widthSingle, heightSingle], 'Parent', h_fig);
ax.ax8 = subplot('Position', [gL0, gb2, widthSingle, heightminus], 'Parent', h_fig);
ax.ax4 = subplot('Position', [gL0, gb1, widthSingle, heightSingle], 'Parent', h_fig);
ax.ax6 = subplot('Position', [gL0, gb0, widthSingle, heightSingle], 'Parent', h_fig);
ax.ax3 = subplot('Position', [gL1, gb3, widthSingle, heightSingle], 'Parent', h_fig);
ax.ax9 = subplot('Position', [gL1, gb2, widthSingle, heightminus], 'Parent', h_fig);
ax.ax5 = subplot('Position', [gL1, gb1, widthSingle, heightSingle], 'Parent', h_fig);
ax.ax7 = subplot('Position', [gL1, gb0, widthSingle, heightSingle], 'Parent', h_fig);

ax.axn1 = subplot('Position', [gL2, gb4, widthDouble, heightplus], 'Parent', h_fig);
ax.axn2 = subplot('Position', [gL2, gb3, widthSingle, heightSingle], 'Parent', h_fig);
ax.axn8 = subplot('Position', [gL2, gb2, widthSingle, heightminus], 'Parent', h_fig);
ax.axn4 = subplot('Position', [gL2, gb1, widthSingle, heightSingle], 'Parent', h_fig);
ax.axn6 = subplot('Position', [gL2, gb0, widthSingle, heightSingle], 'Parent', h_fig);
ax.axn3 = subplot('Position', [gL3, gb3, widthSingle, heightSingle], 'Parent', h_fig);
ax.axn9 = subplot('Position', [gL3, gb2, widthSingle, heightminus], 'Parent', h_fig);
ax.axn5 = subplot('Position', [gL3, gb1, widthSingle, heightSingle], 'Parent', h_fig);
ax.axn7 = subplot('Position', [gL3, gb0, widthSingle, heightSingle], 'Parent', h_fig);

rightBottom = gb0 + 0.04*heightSingle;
rightTop = 1 - buffer;
rightH = rightTop - rightBottom;
rightGap = 0.018;
spH = 0.44 * rightH;
sumH = 0.18 * rightH;
rowH = (rightH - spH - sumH - 3*rightGap) / 2;

spY = rightTop - spH;
sumY = spY - rightGap - sumH;
row1Y = sumY - rightGap - rowH;
row0Y = rightBottom;

ax.axsp = subplot('Position', [gL4, spY, widthDouble, spH], 'Parent', h_fig);
ax.axsum = subplot('Position', [gL4, sumY, widthDouble, sumH], 'Parent', h_fig);
ax.axv1 = subplot('Position', [gL4, row1Y, widthSingle, rowH], 'Parent', h_fig);
ax.axv2 = subplot('Position', [gL5, row1Y, widthSingle, rowH], 'Parent', h_fig);
ax.axv3 = subplot('Position', [gL4, row0Y, widthSingle, rowH], 'Parent', h_fig);
ax.axv4 = subplot('Position', [gL5, row0Y, widthSingle, rowH], 'Parent', h_fig);
end

function updateLabels(src,evt, ax1)
hZ = zoom(src);
hZ.ModeHandle.Blocking = false;
hP = pan(src);
hP.ModeHandle.Blocking = false;
yl = ylim(evt.Axes);
txtObjs = findobj(ax1, 'Type', 'text', 'Tag', 'xregionLabel');
for i = 1:size(txtObjs,1)
    currentPos = txtObjs(i).Position;
    txtObjs(i).Position = [currentPos(1), yl(2) * 0.95, currentPos(3)];
end
hZ.Enable = 'off';
hP.Enable = 'on';
end
