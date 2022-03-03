function [] = PlotBar(x,y,varargin)
%PlotBar Default function to plot bar with individual datapoints
%   input arguments
%       (1) x: location on x-axis
%       (2) y: vector with datapoints

% default properties
Cs = [0.6 0.6 0.6]; % default color
mk = 3; % default marker size

% unpack variable input arguments (optional inputs
if ~isempty(varargin)
    Cs = varargin{1};
    if length(varargin)>1
        mk = varargin{2};
    end
end

%% bar plot
b = bar(x,nanmean(y)); hold on;
b.FaceColor = Cs; b.EdgeColor = Cs;

%% individual datapoitns on top
%   spread datapoint along x-axis
n = length(y);
xrange = 0.2;
dx = (1:n)./n.*xrange - 0.5*xrange + x;
plot(dx,y,'o','MarkerFaceColor',Cs,'Color',[0 0 0],'MarkerSize',mk);





end

