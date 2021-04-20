function configPlot(varargin)
%% configure default plot parameter
% The parameter meanings (and suggested defaults) are::
% 
%   FontSize      = 16   
%   FontName      = 'Times New Roman'
%   Position      = [2, 2, 12, 7.8]   
%   Units         = 'centimeters'
%   LineWidth     = 2
%   AxesLineWidth = 1.2
%
%   Ting Ye
%   07/17 2020


% --- default value
DefaultFontSize = 16;
DefualtAxesLineWidth = 1.2;
DefaultFontName = 'Times New Roman';
DefaultLineWidth = 2;
DefaultFigurePosition = [2, 2, 12, 7.8];
DefaultFigureUnits = 'centimeters';
DefaultMarkerSize = 6;
DefaultColorOrder = getcolor("python");

% parse input parameter
p = inputParser;
addParameter(p, 'FontSize', DefaultFontSize);
addParameter(p, 'FontName', DefaultFontName);
addParameter(p, 'Position', DefaultFigurePosition);
addParameter(p, 'Units', DefaultFigureUnits);
addParameter(p, 'AxesLineWidth', DefualtAxesLineWidth);
addParameter(p, 'LineWidth', DefaultLineWidth);
addParameter(p, 'MarkerSize', DefaultMarkerSize);
addParameter(p, 'ColorOrder', DefaultColorOrder);
parse(p, varargin{:});

% --- figure
set(groot, 'defaultFigureUnits', p.Results.Units);
set(groot, 'defaultFigurePosition', p.Results.Position);
set(groot, 'defaultFigureColor', 'white')

% --- axes (x, y and z)
% set(groot, 'defaultAxesUnits', p.Results.Units);
set(groot, 'defaultAxesFontName', p.Results.FontName);
set(groot, 'defaultAxesFontSize', p.Results.FontSize);
set(groot, 'defaultAxesFontSizeMode', 'manual');
set(groot, 'defaultAxesLineWidth', p.Results.AxesLineWidth);
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1);
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1);
set(groot, 'DefaultAxesColorOrder', p.Results.ColorOrder);

% --- line
set(groot, 'defaultLineLineWidth', p.Results.LineWidth);
set(groot, 'defaultLineMarkerSize', p.Results.MarkerSize);

% --- ConstantLine
set(groot, 'defaultConstantLineLineWidth', p.Results.LineWidth);
set(groot, 'defaultConstantLineFontSize', p.Results.FontSize);

% --- text
set(groot, 'DefaultTextFontSize', p.Results.FontSize);
set(groot, 'DefaultTextFontSizeMode', 'manual');
set(groot, 'DefaultTextFontName', p.Results.FontName);

% --- legend
set(groot, 'defaultLegendFontSize', p.Results.FontSize);
set(groot, 'defaultLegendFontSizeMode', 'manual');
set(groot, 'defaultLegendFontName', p.Results.FontName);
set(groot, 'defaultLegendLocation', 'Best')
set(groot, 'defaultLegendBox', 'off');


end