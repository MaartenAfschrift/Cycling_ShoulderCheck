%% Example on how to analyse data
%--------------------------------

clear all; close all; clc;

% DataPath  = 'S:\Data\fietsproef\Data';
DataPath  = 'E:\fietsproef\Data';
% FigPath = fullfile(DataPath,'Figures\ShoulderCheck');
nPP = 81;
Folders = {'Classic','EBike'};
OrderMeas = {'normal','slow','DualTask','extra','extra2','extra3'};
SensorLocation  = {'Steer','Frame','Trunk','KneeL','KneeR','Pelvis'};

% here we want to compute the variation in steer angle during the "narrow
% lane" task (i.e. first task cycling parcours
ParcousSelected = 'small';

% load the info of the folders
load(fullfile(DataPath,'MatData','ppInfo.mat'),'ppYoung','ppEld');

% flow control
ComputeDataMatrix = true; % Run part to compute datamatrix (or load saved .mat file)
BoolPlot2 = false; % individual plot for each subject with the input data to copumpute variance

% tresholds
threshold_drift = 0.3; % 0.3 radians in this task

figPath = fullfile(pwd,'FigsPaper');


%% select subject
s = 10; % 37 uit zadel , 10 niet uit zadel
f = 1;
i = 1;



%% Get the Datamatrix


ppPath = ['pp_' num2str(s)];
% detect if this is a young or an older subject
BoolEld = any(ppEld ==s);
BoolYoung = any(ppYoung ==s);
% load the data
OutName = [OrderMeas{i} '_data.mat'];
OutPathMat = fullfile(DataPath,'MatData',ppPath,Folders{f});
filename = fullfile(OutPathMat,OutName);
load(filename);
% we don't use the first 3 and last 5 seconds in the
% movement
Rtorso = Phases.Trunk.DualTask.R;
ttorso = Phases.Trunk.DualTask.t;
Rpelvis = Phases.Pelvis.DualTask.R;
tpelvis = Phases.Pelvis.DualTask.t;
Rframe = Phases.Frame.DualTask.R;
tframe = Phases.Frame.DualTask.t;
% get the euler angles
[eulTorso] = GetEulAngles_ShoulderCheck(Rtorso);
[eulpelvis] = GetEulAngles_ShoulderCheck(Rpelvis);
[eulframe] = GetEulAngles_ShoulderCheck(Rframe);
% interpolate eueler angles
eulTorso_int = interp1(ttorso,eulTorso,tframe);
eulPelvis_int = interp1(tpelvis,eulpelvis,tframe);
% relative angles
Q_TorsoFrame = eulTorso_int - eulframe;
Q_PelvisFrame = eulPelvis_int - eulframe;
Q_TorsoPelvis = eulTorso_int -eulPelvis_int;

% einde trial als R frame een bepaalde
% hoek over gaat
% get index turned
t0 = Events.ShoulderCheck(1) - 0.5;
tend = Events.ShoulderCheck(2) + 0.5;
iSel = find(ttorso>t0 & ttorso<tend);
[MinQ,iMin] = min(Q_TorsoFrame(iSel,1));
[MaxQ,iMax] = max(Q_TorsoFrame(iSel,1));
ROM = (MaxQ - MinQ)*180/pi;

[MinQ2,iMin] = min(Q_PelvisFrame(iSel,1));
[MaxQ2,iMax] = max(Q_PelvisFrame(iSel,1));
ROM2 = (MaxQ2 - MinQ2)*180/pi;

[MinQ3,iMin] = min(Q_TorsoPelvis(iSel,1));
[MaxQ3,iMax] = max(Q_TorsoPelvis(iSel,1));
ROM3 = (MaxQ3 - MinQ3)*180/pi;


%% Figure: sensor orientations
figure();

[CPelv, CTors, CFrame] = GetColorsSensorLocation();
lw = 3;


% sensor orientations
plot(tpelvis-t0,eulpelvis(:,1)*180/pi,'Color',CPelv,'LineWidth',lw); hold on;
plot(ttorso-t0,eulTorso(:,1)*180/pi,'Color',CTors,'LineWidth',lw);
plot(tframe-t0,eulframe(:,1)*180/pi,'Color',CFrame,'LineWidth',lw);
set(gcf,'Position',[  979   523   559   274]);

% relative angles
% vline(t0-t0,'k');
% vline(tend-t0,'k');
% title(['subject ' num2str(s) ' ' Folders{f} ' ' OrderMeas{i}]);
% ylabel('orientations');
set(gca,'LineWidth',2);
set(gca,'FontSize',14);
set(gca,'XLim',[0 4.8]);
% set(gca,'YLim',[-1 1.7]);

delete_box
saveas(gcf,fullfile(figPath,'Figure1_SensorOr.svg'),'svg');

%% Figure relative sensor angles

[CTorsoFrame, CTorsoPelvis, CPelvisFrame] = GetColorsJoints();

figure();
plot(tframe-t0,Q_TorsoFrame(:,1)*180/pi,'Color',CTorsoFrame,'LineWidth',lw); hold on;
plot(tframe-t0,Q_PelvisFrame(:,1)*180/pi,'Color',CPelvisFrame,'LineWidth',lw); hold on;
plot(tframe-t0,Q_TorsoPelvis(:,1)*180/pi,'Color',CTorsoPelvis,'LineWidth',lw); hold on;
set(gca,'LineWidth',2);
set(gca,'FontSize',14);
set(gca,'XLim',[0 4.8]);
% set(gca,'YLim',[-1 1.7]);
set(gcf,'Position',[  979   523   559   274]);
delete_box
saveas(gcf,fullfile(figPath,'Figure1_SensorAngle.svg'),'svg');

%% angular velocity sensor
[CPelv, CTors, CFrame] = GetColorsSensorLocation();

figure();
plot(ttorso-t0,Phases.Trunk.DualTask.QdWorld(:,3)*180/pi,'Color',CTors,'LineWidth',lw); hold on;
plot(tframe-t0,Phases.Frame.DualTask.QdWorld(:,3)*180/pi,'Color',CFrame,'LineWidth',lw);
plot(tpelvis-t0,Phases.Pelvis.DualTask.QdWorld(:,3)*180/pi,'Color',CPelv,'LineWidth',lw);
set(gca,'LineWidth',2);
set(gca,'FontSize',14);
set(gca,'XLim',[0 4.8]);
set(gcf,'Position',[  979   523   559   274]);
delete_box
saveas(gcf,fullfile(figPath,'Figure1_qd.svg'),'svg');

%% ROM
figure();
b = bar(1,ROM); b.FaceColor = CTorsoFrame; hold on;
b = bar(2,ROM2); b.FaceColor = CPelvisFrame; hold on;
b = bar(3,ROM3); b.FaceColor = CTorsoPelvis; hold on;
set(gca,'LineWidth',2);
set(gca,'FontSize',14);
set(gca,'XTick',[]);
delete_box
set(gcf,'Position',[547   682   693   296]);
% set(gca,'YLim',[-1 1.7]);
saveas(gcf,fullfile(figPath,'Figure1_ROM.svg'),'svg');

%% steering angle

figure();
isel = find(Phases.SteerAngle.DualTask.t>t0 & Phases.SteerAngle.DualTask.t <tend);
plot(Phases.SteerAngle.DualTask.t-t0,Phases.SteerAngle.DualTask.qSteer(:,1)*180/pi,'Color',CTors,'LineWidth',lw);
set(gca,'XLim',[0 4.8]);
qVar = std(Phases.SteerAngle.DualTask.qSteer(isel,1)*180/pi);

set(gca,'LineWidth',2);
set(gca,'FontSize',14);
delete_box
set(gcf,'Position',[  979   523   559   274]);

saveas(gcf,fullfile(figPath,'Fig1_SteerAngle.svg'),'svg');

disp(['variance in steering angle is ' num2str(qVar)]);
