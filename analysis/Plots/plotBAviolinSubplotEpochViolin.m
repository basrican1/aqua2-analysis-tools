function out = plotBAviolinSubplotEpochViolin(values, t0, datTxt, timeFrameRangeEpochs, varargin)
% Plot BAviolin epoch-grouped violin subplot.
% Can be used standalone or in a parent figure/axes.

p = inputParser;
p.StructExpand = false;
addParameter(p, 'figure', [], @(x) isempty(x) || ishghandle(x, 'figure'));
addParameter(p, 'axes', [], @(x) isempty(x) || ishghandle(x, 'axes'));
addParameter(p, 'style', struct(), @isstruct);
parse(p, varargin{:});
style = mergeBWAPlotStyle(getBWAPlotStyle(), p.Results.style);

[h_fig, h_ax] = resolvePlotAxes('figure', p.Results.figure, 'axes', p.Results.axes, ...
    'name', ['BA_EpochViolin_' char(string(datTxt))]);

epMembIdx = epochMembershipLocal(t0, timeFrameRangeEpochs);
nonMembers = isnan(epMembIdx);
vals = values(~nonMembers);
epIdx = epMembIdx(~nonMembers);
epNames = timeFrameRangeEpochs.epochNames;
if size(epNames,2) == 1
    epNames = epNames';
end

figure(h_fig); axes(h_ax); cla(h_ax);
set(h_ax, 'FontSize', style.fontSizeAxis);
if isempty(vals)
    v = [];
    tCat = categorical([]);
else
    tCat = epNames(epIdx);
    v = violinplot(vals, tCat, 'ShowMedian', true, 'ShowMean', true, 'GroupOrder', epNames);
end
ylabel(h_ax, bwaMetricLabelLatex(datTxt, inferLabelMode(datTxt)), 'Interpreter', 'latex');

out = struct('figure', h_fig, 'axes', h_ax, 'violin', v, 'epMembIdx', epMembIdx, 'nonMembers', nonMembers, 'tCat', tCat);
end

function epMembIdx = epochMembershipLocal(t0, timeFrameRangeEpochs)
epMembIdx = NaN(size(t0));
epFrames = timeFrameRangeEpochs.frameRange;
for epochIndex = 1:size(epFrames,1)
    frameRange = epFrames(epochIndex, 1:2);
    isEpochMember = t0 >= frameRange(1) & t0 <= frameRange(2);
    epMembIdx(isEpochMember) = epochIndex;
end
end

function mode = inferLabelMode(datTxt)
dt = lower(strtrim(char(string(datTxt))));
if startsWith(dt, 'norm ')
    mode = 'norm';
else
    mode = 'raw';
end
end
