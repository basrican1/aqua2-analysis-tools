function [h_fig, h_ax] = resolvePlotAxes(varargin)
% Resolve plotting target axes.
% If figure/axis are not provided, creates a new figure with one axes.
%
% Name-value:
%   'figure' : figure handle (optional)
%   'axes'   : axes handle (optional)
%   'name'   : figure name when creating new (default: '')

p = inputParser;
p.StructExpand = false;
addParameter(p, 'figure', [], @(x) isempty(x) || ishghandle(x, 'figure'));
addParameter(p, 'axes', [], @(x) isempty(x) || ishghandle(x, 'axes'));
addParameter(p, 'name', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

h_fig = p.Results.figure;
h_ax = p.Results.axes;

if isempty(h_fig) && isempty(h_ax)
    if strlength(string(p.Results.name)) > 0
        h_fig = figure('Name', char(string(p.Results.name)));
    else
        h_fig = figure();
    end
    h_ax = axes('Parent', h_fig);
    return
end

if isempty(h_ax) && ~isempty(h_fig)
    figure(h_fig);
    h_ax = axes('Parent', h_fig);
    return
end

if ~isempty(h_ax) && isempty(h_fig)
    h_fig = ancestor(h_ax, 'figure');
end
end
