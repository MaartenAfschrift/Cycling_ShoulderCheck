
%% Compute rotaton axis of the steer for each subject
clear all; close all; clc;
% Path information
DataPath  = 'E:\Data\Fietsproef\MatData';
FigPath = 'E:\Data\Fietsproef\Figures\HingeSteer';
nPP = 81;
Folders = {'Classic','EBike'};
OrderMeas = {'normal','slow','DualTask','extra','extra2','extra3'};

% Flow control
BoolGetAxis = false;
BoolGetAngle =  true;

%% Get the axis of rotaton of the steer for each subject
if BoolGetAxis
    %ct = 1;
    BoolPlot = true;
    % loop over all subjects
    for s = 1:nPP
        ppPath = ['pp_' num2str(s)];
        for f = 1:length(Folders)
            i = 1; % callibration is always in file 1, expect for:
            if (s==1) || (s==3 && f==2) % we know that the callibration is in file 2 here
                i = 2;
            end
            % load the data
            OutName = [OrderMeas{i} '_data.mat'];
            OutPathMat = fullfile(DataPath,ppPath,Folders{f});
            filename = fullfile(OutPathMat,OutName);
            if exist(filename,'file')
                % load the mat file (processed with ExampleBatch2)
                load(filename,'tTrigger','Data','header');
                % compute the orientation of the steer axis in the files
                % with callibration motion
                % test if this file is a callibration file
                if length(tTrigger) == 12
                    
                    % get the phase with the callibration procedure
                    dt = [tTrigger(2) tTrigger(3)];
                    timeSpan = getCallibrationPhase(Data.Steer.t,dt,Data.Steer.QdWorld);
                    % we want to start in "neutral position". THis is 2seconds
                    % before start callibration motion
                    timeSpan(1) = timeSpan(1)-3;
                    % interpolate rotation matrices if needed
                    if (length(Data.Frame.t) ~= length(Data.Steer.t)) || (any((Data.Frame.t-Data.Steer.t)~=0))
                        [Data.Frame.Rint, Data.Steer.Rint, tint] = InterpolateRotMatrices(Data.Frame.R,...
                            Data.Steer.R,Data.Frame.t,Data.Steer.t);
                        disp(['Interpolated rotation matrices for file: ' filename]);
                    else
                        Data.Frame.Rint = Data.Frame.R;
                        Data.Steer.Rint = Data.Steer.R;
                        tint = Data.Frame.t;
                    end
                    % get hinge
                    [Rax,n_steer,n_frame] = GetHingeAxis(Data.Frame.Rint,Data.Steer.Rint,tint,timeSpan,BoolPlot);
                    save(fullfile(OutPathMat,'RotAxis_Steer.mat'),'Rax','n_steer','n_frame');
                    if BoolPlot
                        OutPathFig = fullfile(FigPath,ppPath,Folders{f});
                        if ~isfolder(OutPathFig)
                            mkdir(OutPathFig);
                        end
                        saveas(gcf,fullfile(OutPathFig,[OrderMeas{i} '.fig']));
                        close(gcf);
                    end
                else
                    disp(['Could not find callibration trial in subject: ' num2str(s) ' - ' Folders{f}]);
                end
            end
        end
    end
end

%% Compute the rotation in each file
if BoolGetAngle
    % loop over all subjects
    for s = 1:nPP
        ppPath = ['pp_' num2str(s)];
        for f = 1:length(Folders)
            OutPathMat = fullfile(DataPath,ppPath,Folders{f});
            fileSteer = fullfile(OutPathMat,'RotAxis_Steer.mat');
            for i =1:3
                % load the data
                OutName = [OrderMeas{i} '_data.mat'];
                filename = fullfile(OutPathMat,OutName);
                if exist(filename,'file') && exist(fileSteer,'file')
                    
                    % load the datafiles
                    load(filename,'tTrigger','Data','header');
                    load(fileSteer,'Rax');
                    
                    % compute the steer angle
                    if ~isempty(Data)
                        if  (length(Data.Frame.t) ~= length(Data.Steer.t)) || (any((Data.Frame.t-Data.Steer.t)~=0))
                            [Data.Frame.Rint, Data.Steer.Rint, tint] = InterpolateRotMatrices(Data.Frame.R,Data.Steer.R,Data.Frame.t,Data.Steer.t);
                            disp(['Interpolated rotation matrices for file: ' filename]);
                        else
                            Data.Frame.Rint = Data.Frame.R;
                            Data.Steer.Rint = Data.Steer.R;
                            tint = Data.Frame.t;
                        end
                        [SteerAngle.q] = GetAngleSteer(Data.Frame.Rint,Data.Steer.Rint,Rax);
                        SteerAngle.t = tint;
                        
                        % append the matlab structure with this new information
                        save(filename,'SteerAngle','-append');
                    end
                end
            end
        end
        disp(['Subject ' num2str(s) ' / ' num2str(nPP)]);
    end
end






