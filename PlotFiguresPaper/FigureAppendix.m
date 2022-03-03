%% Bar Plots ROM and Steering angle

% Figure with IMU outcomes during cycling parcours
% clean the matlab environment
clear all; close all; clc;

% path information
[TempPath,~,~] = fileparts(mfilename('fullpath'));
MainPath = TempPath(1:end-17);
SensorOr = load(fullfile(MainPath,'Outcomes','ShouldCheckROM.mat'),'DataMatrix','header_DataMatrix');

% path to save figures
figPath = fullfile(MainPath,'Outcomes','FigsPaper');

%% FIgure steering angle

h = figure();
set(h,'Position',[113   394   765   623]);

CYoung = [0 114 178]./255;
CEld = [213 94 0]./255;
mk = 3;

% plot figure
iBike = [1 1 2 2];
ispeed = [1 2 1 2];

for i=1:4
    DataMatrix = SensorOr.DataMatrix;
    subplot(2,2,1);
    iSelY = DataMatrix(:,5) == 0 &  DataMatrix(:,3) == ispeed(i) & DataMatrix(:,2) == iBike(i) & DataMatrix(:,6) == 0;
    PlotBar(i*3-2,DataMatrix(iSelY,9),CYoung,mk); hold on;
    iSelE = DataMatrix(:,5) == 1 &  DataMatrix(:,3) == ispeed(i) & DataMatrix(:,2) == iBike(i) & DataMatrix(:,6) == 0;
    PlotBar(i*3-1,DataMatrix(iSelE,9),CEld,mk); hold on;  
    set(gca,'Box','off');
    title('steering angle')
end
ylabel('standard deviation [deg]')


DataMatrix = SensorOr.DataMatrix;
iCol = [4 7 8];
TitleSel = {'ROM Frame-C7','ROM Frame-Pelvis','ROM Torso-Pelvis'};

for i=1:3
    for j =1:4
        subplot(2,2,i+1);
        iSelY = DataMatrix(:,5) == 0 &  DataMatrix(:,3) == ispeed(j) & DataMatrix(:,2) == iBike(j) & DataMatrix(:,6) == 0;
        PlotBar(j*3-2,DataMatrix(iSelY,iCol(i)),CYoung,mk); hold on;
        iSelE = DataMatrix(:,5) == 1 &  DataMatrix(:,3) == ispeed(j) & DataMatrix(:,2) == iBike(j) & DataMatrix(:,6) == 0;
        PlotBar(j*3-1,DataMatrix(iSelE,iCol(i)),CEld,mk); hold on;
    end    
    ylabel('ROM [deg]');    
    set(gca,'Box','off');
    title(TitleSel{i})
end

for i=1:4
    subplot(2,2,i)
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1.5);
    set(gca,'XTick',[1.5 4.5 7.5 10.5])
    set(gca,'XTickLabel',{'N-Bike','N-Bike Slow','E-Bike','E-Bike Slow'})
    set(gca,'XTickLabelRotation',60)
end

