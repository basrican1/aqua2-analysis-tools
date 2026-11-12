function txt = bwaMetricLabelLatex(datTxt, mode)
% Build consistent LaTeX labels for BAviolin metric axes.
% mode: 'raw', 'norm', 'cell_raw', 'cell_norm'

if nargin < 2 || isempty(mode)
    mode = 'raw';
end
mode = lower(char(string(mode)));

raw = char(string(datTxt));
raw = strtrim(raw);
if startsWith(lower(raw), 'norm ')
    raw = strtrim(raw(6:end));
end
if strlength(string(raw)) >= 2 && strcmp(raw(end-1:end), '_c')
    raw = raw(1:end-2);
end

core = metricCore(raw); % latex fragment WITHOUT outer $...$
switch mode
    case 'raw'
        txt = ['$' core '$'];
    case 'norm'
        txt = ['$\log_2\mathrm{FC}\left(' core '\right)$'];
    case 'cell_raw'
        txt = ['$' core '_{\mathrm{cell}}$'];
    case 'cell_norm'
        txt = '$\log_2\mathrm{FC}_{\mathrm{cell}}$';
    otherwise
        txt = ['$' core '$'];
end
end

function core = metricCore(raw)
r = lower(strrep(strrep(strrep(raw, '_', ''), '.', ''), ' ', ''));
switch r
    case {'dffauc','curvedffauc','dffaucsec'}
        core = '\mathrm{dffAUC}\ (\Delta F/F\cdot \mathrm{s})';
    case {'dffmax2','curvedffmax2'}
        core = '\mathrm{dffMax2}';
    case {'durationframes','durationsec','curveduration'}
        core = '\mathrm{Duration}\ (\mathrm{s})';
    case {'rise19frames','rise19sec','curverise19'}
        core = '\mathrm{Rise}_{19}\ (\mathrm{s})';
    case {'fall91frames','fall91sec','curvefall91'}
        core = '\mathrm{Fall}_{91}\ (\mathrm{s})';
    case {'decaytauframes','decaytausec','curvedecaytau'}
        core = '\tau_{\mathrm{decay}}\ (\mathrm{s})';
    case {'areapx2','areaum2','basicarea'}
        core = '\mathrm{Area}\ (\mu\mathrm{m}^2)';
    case {'perimeterpx','perimeterum'}
        core = '\mathrm{Perimeter}\ (\mu\mathrm{m})';
    case {'propagationdistancepx','propagationdistanceum','propagationpropdist'}
        core = '\mathrm{Propagation\ Distance}\ (\mu\mathrm{m})';
    case {'meanneighborentrydelayframes','meanneighborentrydelaysec'}
        core = '\mathrm{Neighbor\ Entry\ Delay}\ (\mathrm{s})';
    case {'meansomacentroiddistpx','meansomacentroiddistum'}
        core = '\mathrm{Soma\ Centroid\ Distance}\ (\mu\mathrm{m})';
    case {'meansomastartdistpx','meansomastartdistum'}
        core = '\mathrm{Soma\ Start\ Distance}\ (\mu\mathrm{m})';
    case {'meansomaenddistpx','meansomaenddistum'}
        core = '\mathrm{Soma\ End\ Distance}\ (\mu\mathrm{m})';
    case {'meanpropagationdistancepx','meanpropagationdistanceum'}
        core = '\mathrm{Mean\ Propagation\ Distance}\ (\mu\mathrm{m})';
    case {'meanperimeterpx','meanperimeterum'}
        core = '\mathrm{Mean\ Perimeter}\ (\mu\mathrm{m})';
    case {'maxpropspeed','propagationmaxpropspeed'}
        core = '\mathrm{Max\ Prop\ Speed}';
    otherwise
        safe = strrep(raw, '_', '\_');
        core = ['\mathrm{' safe '}'];
end
end
