function files = plotPost2FrequencyFigures(result, figureSpecPath, outputDir)
%PLOTPOST2FREQUENCYFIGURES Render factor-aware canonical frequency figures.
% Figures are descriptive. Connected lines represent repeated epochs within a
% cell or recording; colors and requested comparisons come from figure_spec.

if ~(isstruct(result) && isfield(result, 'canonical') && isfield(result.canonical, 'frequency') && ...
        isfield(result.canonical.frequency, 'ratesByCellEpoch'))
    error('plotPost2FrequencyFigures:MissingFrequency', 'Result has no canonical cell-by-epoch frequency table.');
end
T = result.canonical.frequency.ratesByCellEpoch;
if isempty(T), error('plotPost2FrequencyFigures:EmptyFrequency', 'Frequency table is empty.'); end
spec = jsondecode(fileread(char(figureSpecPath)));
style = resolveStyle(spec);
if exist(outputDir, 'dir') ~= 7, mkdir(outputDir); end

files = strings(0,1);
panels = structList(spec.comparison_panels);
f = figure('Visible','off','Color','w','Position',[100 80 1280 390*numel(panels)]);
tiledlayout(numel(panels),1,'TileSpacing','compact','Padding','compact');
axesToShare = gobjects(numel(panels),1);
for pp = 1:numel(panels)
    ax = nexttile; axesToShare(pp) = ax; hold(ax,'on');
    P = panels(pp); groups = comparisonGroups(P);
    handles = gobjects(numel(groups),1);
    for gg = 1:numel(groups)
        G = groups(gg);
        rows = string(T.bath_condition_id) == G.bath_condition_id & string(T.stim_pattern) == G.stim_pattern;
        if ~any(rows), error('plotPost2FrequencyFigures:MissingGroup', 'No cell rows for %s.', G.label); end
        color = groupColor(style, G.bath_condition_id, G.stim_pattern);
        [x, Y] = cellEpochMatrix(T(rows,:));
        faded = 0.62 * color + 0.38 * [1 1 1];
        plot(ax, x, Y', '-', 'Color', faded, 'LineWidth', 0.7);
        mu = mean(Y,1,'omitnan'); n = sum(isfinite(Y),1); sem = std(Y,0,1,'omitnan') ./ sqrt(max(n,1));
        good = isfinite(mu) & isfinite(sem);
        if any(good)
            xx = x(good); fill(ax,[xx fliplr(xx)],[mu(good)+sem(good) fliplr(mu(good)-sem(good))], ...
                color,'FaceAlpha',style.fill_alpha,'EdgeColor','none');
        end
        handles(gg) = plot(ax,x,mu,'-o','Color',color,'MarkerFaceColor',color, ...
            'LineWidth',style.line_width,'MarkerSize',6);
    end
    yline(ax,0,':','Color',[.45 .45 .45]);
    xlim(ax,[.7 3.3]); xticks(ax,1:3); xticklabels(ax,epochLabels(T));
    ylabel(ax,'Events per minute per cell'); title(ax,char(string(P.title)),'Interpreter','none');
    legend(ax,handles,string({groups.label}),'Location','best','Interpreter','none'); box(ax,'off');
end
applySharedYLimits(axesToShare, style); sgtitle('Canonical post-2 cell frequency: mean ± SEM with within-cell epoch trajectories');
files(end+1,1) = exportFigure(f, outputDir, 'post2_frequency_requested_comparisons'); close(f);

f = figure('Visible','off','Color','w','Position',[100 100 1350 650]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
R = recordingSummary(T);
groups = unique(string(R.bath_condition_id) + " | " + string(R.stim_pattern),'stable');
ax = nexttile; hold(ax,'on'); handles = gobjects(numel(groups),1);
for gg = 1:numel(groups)
    rows = string(R.bath_condition_id) + " | " + string(R.stim_pattern) == groups(gg);
    color = groupColor(style, R.bath_condition_id(find(rows,1)), R.stim_pattern(find(rows,1)));
    [x,Y] = recordingEpochMatrix(R(rows,:));
    plot(ax,x,Y','-','Color',0.62*color+0.38,'LineWidth',0.9);
    handles(gg) = plot(ax,x,mean(Y,1,'omitnan'),'-o','Color',color,'MarkerFaceColor',color,'LineWidth',2.3);
end
yline(ax,0,':','Color',[.45 .45 .45]); xlim(ax,[.7 3.3]); xticks(ax,1:3); xticklabels(ax,epochLabels(T));
ylabel(ax,'Mean events per minute per cell'); title(ax,'Recording-level trajectories'); legend(ax,handles,replace(groups,' | ',' '),'Location','best','Interpreter','none'); box(ax,'off');

ax = nexttile; hold(ax,'on'); E = R(R.epoch_name=="post_early",:); L = R(R.epoch_name=="post_late",:);
[matched, idx] = ismember(string(E.recording_id),string(L.recording_id)); E=E(matched,:); L=L(idx(matched),:);
allVals = [E.mean_rate_events_per_min;L.mean_rate_events_per_min]; lim = paddedLimits(allVals);
for gg = 1:numel(groups)
    rows = string(E.bath_condition_id) + " | " + string(E.stim_pattern) == groups(gg);
    color = groupColor(style,E.bath_condition_id(find(rows,1)),E.stim_pattern(find(rows,1)));
    scatter(ax,E.mean_rate_events_per_min(rows),L.mean_rate_events_per_min(rows),52,'filled', ...
        'MarkerFaceColor',color,'MarkerEdgeColor','w');
end
plot(ax,lim,lim,':k'); xline(ax,0,':','Color',[.45 .45 .45]); yline(ax,0,':','Color',[.45 .45 .45]);
xlim(ax,lim); ylim(ax,lim); axis(ax,'square'); xlabel(ax,'Post early mean frequency'); ylabel(ax,'Post late mean frequency'); title(ax,'Recording-level early versus late'); box(ax,'off');
sgtitle('Canonical post-2 recording summaries');
files(end+1,1) = exportFigure(f, outputDir, 'post2_frequency_recording_relationships'); close(f);
end

function style = resolveStyle(spec)
root = char(string(spec.shared_metadata_root));
adapter = fullfile(root,'adapters','matlab');
if exist(adapter,'dir') ~= 7, error('plotPost2FrequencyFigures:SharedMetadata','Missing shared palette adapter: %s',adapter); end
addpath(adapter); [rgb,palette] = getSharedPalette(char(string(spec.palette))); %#ok<ASGLU>
style = struct('rgb',rgb,'bath',string(spec.factor_levels.bath_condition_id(:)), ...
    'stim',string(spec.factor_levels.stim_pattern(:)), ...
    'fill_alpha',double(spec.display.summary_fill_alpha), ...
    'line_width',double(spec.display.summary_line_width), ...
    'share_y',strcmpi(string(spec.axes.share_y_limits),'figure'), ...
    'padding',double(spec.axes.y_limit_padding_fraction));
end

function groups = comparisonGroups(panel)
levels = string(panel.level_order(:)); labels = string(panel.group_labels(:)); varying = string(panel.varying_factor); fixed = panel.fixed_factors;
groups = repmat(struct('bath_condition_id',"",'stim_pattern',"",'label',""),numel(levels),1);
for ii=1:numel(levels)
    if varying=="bath_condition_id", bath=levels(ii); stim=string(fixed.stim_pattern); else, bath=string(fixed.bath_condition_id); stim=levels(ii); end
    groups(ii)=struct('bath_condition_id',bath,'stim_pattern',stim,'label',labels(ii));
end
end

function color = groupColor(style,bath,stim)
bi=find(style.bath==string(bath),1); si=find(style.stim==string(stim),1);
if isempty(bi)||isempty(si), error('plotPost2FrequencyFigures:ColorLevel','Undeclared color level.'); end
color=style.rgb((bi-1)*numel(style.stim)+si,:);
end

function [x,Y] = cellEpochMatrix(T)
x = 1:3; ids=unique(string(T.cell_uid),'stable'); Y=nan(numel(ids),3);
for ii=1:numel(ids), rows=string(T.cell_uid)==ids(ii); for ee=1:3, v=T.rate_events_per_min(rows & double(T.epoch_idx)==ee); if isscalar(v),Y(ii,ee)=v;end,end,end
end

function R = recordingSummary(T)
keys=unique(T(:,{'recording_id','bath_condition_id','stim_pattern','epoch_idx','epoch_name'}),'rows','stable');
R=keys; R.mean_rate_events_per_min=nan(height(R),1);
for ii=1:height(R), rows=string(T.recording_id)==string(R.recording_id(ii)) & double(T.epoch_idx)==double(R.epoch_idx(ii)); R.mean_rate_events_per_min(ii)=mean(T.rate_events_per_min(rows),'omitnan'); end
end

function [x,Y] = recordingEpochMatrix(R)
x=1:3; ids=unique(string(R.recording_id),'stable'); Y=nan(numel(ids),3);
for ii=1:numel(ids), rows=string(R.recording_id)==ids(ii); for ee=1:3,v=R.mean_rate_events_per_min(rows & double(R.epoch_idx)==ee);if isscalar(v),Y(ii,ee)=v;end,end,end
end

function labels = epochLabels(T)
labels = strings(1,3); for ee=1:3, name=string(T.epoch_name(find(double(T.epoch_idx)==ee,1))); labels(ee)=replace(name,{'pre_stim_baseline','post_early','post_late'},{'Pre-stim baseline','Post early','Post late'}); end
end
function applySharedYLimits(axs,style), if ~style.share_y,return,end; axs=axs(isgraphics(axs)); lo=min(arrayfun(@(a)a.YLim(1),axs)); hi=max(arrayfun(@(a)a.YLim(2),axs)); pad=style.padding*(hi-lo);if pad==0,pad=1;end;set(axs,'YLim',[lo-pad hi+pad]);end
function lim=paddedLimits(x),x=x(isfinite(x));if isempty(x),lim=[-1 1];return,end;pad=.08*range(x);if pad==0,pad=1;end;lim=[min(x)-pad max(x)+pad];end
function file=exportFigure(f,out,stem),file=fullfile(out,stem+".png");exportgraphics(f,file,'Resolution',180);savefig(f,fullfile(out,stem+".fig"));end
function out=structList(value),if ~iscell(value),out=value;return,end;names=strings(0,1);for ii=1:numel(value),names=unique([names;string(fieldnames(value{ii}))],'stable');end;out=repmat(cell2struct(cell(numel(names),1),cellstr(names),1),numel(value),1);for ii=1:numel(value),f=fieldnames(value{ii});for jj=1:numel(f),out(ii).(f{jj})=value{ii}.(f{jj});end,end,end
