function style = getBWAPlotStyle()
% Central plotting style defaults for BWA plotting functions.
% Pass a partial override struct through 'style' name-value where supported.

style = struct();
style.figurePositionBA = [100 100 1600 800];
style.figurePositionPolars = [50 50 1200 400];
style.fontSizeAxis = 10;
style.fontSizeEpochLabel = 10;
style.fontSizeLegend = 8;
style.legendLocation = 'southoutside';
style.tileSpacing = 'none';
style.tilePadding = 'compact';

style.lineWidthMain = 1.5;
style.lineWidthBaseline = 2.5;
style.lineWidthError = 1.0;
style.markerSize = 10;
style.scatterSize = 8;

style.colors = struct();
style.colors.groupA = [0 0 1];
style.colors.groupC = [1 0 0];
style.colors.baseline = [0 0 0];
style.colors.text = [0 0 0];
end
