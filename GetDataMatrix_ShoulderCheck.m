%% Analyzing data shoulder check
%-------------------------------

% start with clean matlab
clear all; close all; clc;

% path information
[MainPath,~] = fileparts(mfilename('fullpath'));
DataPath  = fullfile(MainPath,'Data');

% add functions to the matlab path
addpath(genpath(fullfile(MainPath,'Functions')));

% subject information
nPP = 81;
Folders = {'Classic','EBike'};
OrderMeas = {'normal','Slow'};
SensorLocation  = {'Steer','Frame','Trunk','Pelvis'};

% load the info of the folders
load(fullfile(DataPath,'ppInfo.mat'),'ppYoung','ppEld');

%% Read the excel table with information on the task
FileTaskInfo = fullfile(DataPath,'TaskInformation.xlsx');
[ShoulderCheckInfo] = GetShoulderCheckInfo(FileTaskInfo);
ShoulderCheckInfo.pp = [ShoulderCheckInfo.ppIDYoung; ShoulderCheckInfo.ppIDOld];
ShoulderCheckInfo.data = [ShoulderCheckInfo.DatYoung; ShoulderCheckInfo.DatOlder];
ShoulderCheckInfo.header = ShoulderCheckInfo.HeadersOlder;


%% Get the Datamatrix
DataMatrix = nan(100*3*2,9); % pre allocate matrix with all the data
header_DataMatrix =  {'s-ID','bike-ID','Speed-ID','ROM-FrameTorso','BoolElderly',...
    'Error','ROM-FramePelvis','ROM-PelvisTorso','SteeringAngle','CorrSteerTorso'};
diary('LogExample_ShoulderCheck.txt');
ct = 1;
for s = 1:nPP
    ppPath = ['pp_' num2str(s)];
    % detect if this is a young or an older subject
    BoolEld = any(ppEld ==s);
    BoolYoung = any(ppYoung ==s);
    if BoolEld && BoolYoung
        disp(['Subject ' num2str(s) ' is both young and old. adapt this in the excel file']);
    end
    for f = 1:length(Folders)
        % load the axis of rotation to compute the steering angle (this
        % axis is computed in the script GetRotationAxisSteer_Subjects.m)
        OutPathMat = fullfile(DataPath,ppPath,Folders{f});
        RotAx = load(fullfile(OutPathMat,'RotAxis_Steer.mat'),'Rax','n_steer','n_frame');
        for i =1:length(OrderMeas)
            % load the data
            OutName = [OrderMeas{i} '_data.mat'];            
            filename = fullfile(OutPathMat,OutName);
            if exist(filename,'file')
                load(filename,'Data','GUIvar','Events');
                if ~exist('GUIvar','var')
                    GUIvar = [];
                end
                if ~exist('Events','var')
                    Events = [];
                end
                if exist('Data','var')
                    BoolErrorFlag = 0;
                    % check if task was performed according to
                    % instructions
                    iData = find(ShoulderCheckInfo.pp== s);
                    if strcmp(Folders{f},'Classic')
                        if strcmp(OrderMeas{i},'normal')
                            headerSel = {'ErrorCone-normal-classic',...
                                'FootOnGround-normal-classic'};
                        elseif strcmp(OrderMeas{i},'slow')
                            headerSel = {'ErrorCone-slow-classic',...
                                'FootOnGround-slow-classic'};
                        end
                    elseif strcmp(Folders{f},'EBike')
                        if strcmp(OrderMeas{i},'normal')
                            headerSel = {'ErrorCone-normal-ebike',...
                                'FootOnGround-normal-ebike'};
                        elseif strcmp(OrderMeas{i},'slow')
                            headerSel = {'ErrorCone-slow-ebike',...
                                'FootOnGround-slow-ebike'};
                        end
                    end
                    % select colIndices
                    IndsColInfo = nan(length(headerSel),1);
                    for ic=1:length(headerSel)
                        IndsColInfo(ic) = find(strcmp(ShoulderCheckInfo.header,headerSel{ic}));
                    end
                    ErrorTask = ShoulderCheckInfo.data(iData,IndsColInfo);
                    if sum(ErrorTask) == 0 && isfield(Data,'SteerAngle') && ~isempty(Data.SteerAngle.t)
                        % get the sensor orientations
                        Rtorso = Data.Trunk.R;
                        ttorso = Data.Trunk.t;
                        Rframe = Data.Frame.R;
                        tframe = Data.Frame.t;
                        % we had some issues with the pelvis sensor
                        if ~isfield(Data,'Pelvis')
                            Rpelvis = [];
                            tpelvis = [];
                        else
                            Rpelvis = Data.Pelvis.R;
                            tpelvis = Data.Pelvis.t;
                        end
                        % log errors witht the pelvis sensor
                        BoolPelvisError = false;
                        if isempty(Rpelvis)
                            BoolPelvisError = true;
                        end
                        % The structure Events contains information on the
                        % start and end of the shoulder check (wich was
                        % manualy detected using a GUI (GUI_ShoulderCheckEvents)
                        if (isfield(Events,'ShoulderCheck') && ~isempty(Events.ShoulderCheck) &&   ~any(isnan(Events.ShoulderCheck)))
                            % get the euler angles
                            [eulTorso] = GetEulAngles_ShoulderCheck(Rtorso);
                            if ~BoolPelvisError
                                [eulpelvis] = GetEulAngles_ShoulderCheck(Rpelvis);
                            end
                            [eulframe] = GetEulAngles_ShoulderCheck(Rframe);
                            % interpolate eueler angles
                            eulTorso_int = interp1(ttorso,eulTorso,tframe);
                            if ~BoolPelvisError
                                eulPelvis_int = interp1(tpelvis,eulpelvis,tframe);
                            end
                            % relative angles
                            Q_TorsoFrame = eulTorso_int - eulframe;
                            if ~BoolPelvisError
                                Q_PelvisFrame = eulPelvis_int - eulframe;
                                Q_TorsoPelvis = eulTorso_int -eulPelvis_int;
                            end
                            % get data in selected time window
                            t0 = Events.ShoulderCheck(1) - 0.5; % start 0.5s before detected start
                            tend = Events.ShoulderCheck(2) + 0.5; % end 0.5s after detected end
                            iSel = find(ttorso>t0 & ttorso<tend);
                            [MinQ,iMin] = min(Q_TorsoFrame(iSel,1));
                            [MaxQ,iMax] = max(Q_TorsoFrame(iSel,1));
                            ROM = (MaxQ - MinQ)*180/pi; % convert to deg
                            if isempty(ROM)
                                ROM = NaN;
                            end
                            % compute steering angle from axis
                            if  (length(Data.Frame.t) ~= length(Data.Steer.t)) || (any((Data.Frame.t-Data.Steer.t)~=0))
                                [Data.Frame.Rint, Data.Steer.Rint, tint] = InterpolateRotMatrices(Data.Frame.R,Data.Steer.R,Data.Frame.t,Data.Steer.t);
                                disp(['Interpolated rotation matrices for file: ' filename]);
                            else
                                Data.Frame.Rint = Data.Frame.R;
                                Data.Steer.Rint = Data.Steer.R;
                                tint = Data.Frame.t;
                            end
                            [SteerAngle.q] = GetAngleSteer(Data.Frame.Rint,Data.Steer.Rint,RotAx.Rax);
                            SteerAngle.t = tint;
                            % get standard deviation in steering angle
                            iSelSteer = find(tint>t0 & tint<tend);
                            q = Data.SteerAngle.qSteer(:,1);
                            qVarDeg = std(q(iSelSteer)*180/pi);
                            if qVarDeg >20 % variances above 20 deg are impossible during this task
                                disp(['possible error in file: ' filename ' remove this file from the analysis']);
                                qVarDeg = NaN;
                            end
                            % correlation between steering angle
                            % and orientation of the torso
                            eulTorso_intSteer = interp1(ttorso,eulTorso(:,1),Data.SteerAngle.t(iSelSteer))';
                            qSteer = q(iSelSteer,1);
                            rho = corr(eulTorso_intSteer,qSteer);
                            if isfield(GUIvar,'Shoulder_Drift') && ~GUIvar.Shoulder_Drift && ~BoolPelvisError
                                [MinQ2,iMin] = nanmin(Q_PelvisFrame(iSel,1));
                                [MaxQ2,iMax] = nanmax(Q_PelvisFrame(iSel,1));
                                ROM2 = (MaxQ2 - MinQ2)*180/pi;
                                [MinQ3,iMin] = min(Q_TorsoPelvis(iSel,1));
                                [MaxQ3,iMax] = max(Q_TorsoPelvis(iSel,1));
                                ROM3 = (MaxQ3 - MinQ3)*180/pi;
                            else
                                ROM2 = NaN;
                                ROM3 = NaN;
                            end
                        else
                            ROM = NaN;
                            ROM2 = NaN;
                            ROM3 = NaN;
                            qVarDeg = NaN;
                            rho = NaN;
                        end
                        % store the ROM in the datamatrix
                        DataMatrix(ct,1) = s;
                        DataMatrix(ct,2) = f;
                        DataMatrix(ct,3) = i;
                        DataMatrix(ct,4) =  ROM ;
                        DataMatrix(ct,5) = BoolEld;
                        DataMatrix(ct,6) = BoolErrorFlag;
                        DataMatrix(ct,7) = ROM2;
                        DataMatrix(ct,8) = ROM3;
                        DataMatrix(ct,9) = qVarDeg;
                        DataMatrix(ct,10) = rho;
                        ct = ct+1;
                    end
                end
                clear Data GUIvar Events
            end
        end
    end
    disp(['Subject ' num2str(s) ' / ' num2str(nPP)]);
end
DataMatrix(ct:end,:) = [];
% save the datamatrix
if ~isfolder(fullfile(MainPath,'Outcomes'))
    mkdir(fullfile(MainPath,'Outcomes'));
end
save(fullfile(MainPath,'Outcomes','ShouldCheckROM.mat'),'DataMatrix','header_DataMatrix');
diary off
