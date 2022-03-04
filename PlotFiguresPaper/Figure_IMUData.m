% ----------------------------------------------
% Figure with IMU outcomes during shoulder check
%-----------------------------------------------

% clean the matlab environment
clear all; close all; clc;

% path information
[TempPath,~,~] = fileparts(mfilename('fullpath'));
MainPath = TempPath(1:end-17);
Dat = load(fullfile(MainPath,'Outcomes','ShouldCheckROM.mat'),'DataMatrix','header_DataMatrix');

% location to save the figures
figPath = fullfile(MainPath,'Outcomes','FigsPaper');

% open a diary to log messages
diary('Figure2_Log.txt');    

% open a figure
h = figure();
set(h,'Position',[113   394   765   623]);

% color codes for young, older subjects
CYoung = [0 114 178]./255;
CEld = [213 94 0]./255;
mk = 3;

% get indices for specific columns
iAge = strcmp(Dat.header_DataMatrix,'BoolElderly');
iSpeed = strcmp(Dat.header_DataMatrix,'Speed-ID');
iBike = strcmp(Dat.header_DataMatrix,'bike-ID');
iError = strcmp(Dat.header_DataMatrix,'Error');

% selected data for plots
DataName = {'SteeringAngle','ROM-FrameTorso','ROM-FramePelvis','ROM-PelvisTorso'};

% title for selected data
TitleSel = {'Steering angle','ROM Frame-Torso','ROM Frame-Pelvis','ROM Pelvis-Torso'};

% direct reference to data matrix
DataMatrix = Dat.DataMatrix;

% plot sensor orientations
for i=1:4
    subplot(2,2,i);

    % select col
    iCol = strcmp(Dat.header_DataMatrix,DataName{i});
    
    % plot datapoints
    iSelY = DataMatrix(:,iAge) == 0 & DataMatrix(:,iSpeed) == 1 & ...
        DataMatrix(:,iBike) == 1 & DataMatrix(:,iError) == 0;
    PlotBar(1,DataMatrix(iSelY,iCol),CYoung,mk); hold on;
    iSelE = DataMatrix(:,iAge) == 1 &  DataMatrix(:,iSpeed) == 1 &...
        DataMatrix(:,iBike) == 1 & DataMatrix(:,iError) == 0;
    PlotBar(2,DataMatrix(iSelE,iCol),CEld,mk); hold on;
        
    % paired t-test on the selected data
    [pairedttest,p,ci,stats] = ttest2(DataMatrix(iSelY,iCol),...
        DataMatrix(iSelE,iCol),0.05);

    % display statistics
    disp(TitleSel{i});
    disp(['number of young subjects ' , num2str(sum(~isnan(DataMatrix(iSelY,iCol))))]);
    disp(['number of older subjects ' , num2str(sum(~isnan(DataMatrix(iSelE,iCol))))]);
    disp(' ');
    title([TitleSel{i} ': p = ' num2str(p)]);    
    disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ' , p = ' num2str(p)])
    disp(' ');

    % labels
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1.5);
    if i == 1
        ylabel('Angle [deg]');
    else
        ylabel('ROM [deg]');
    end
    set(gca,'XTick',1:2);
    set(gca,'XTickLabel',{'Young','Older'});
    set(gca,'Box','off')    
end

