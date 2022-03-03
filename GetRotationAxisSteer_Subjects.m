
%% Compute rotaton axis of the steer for each subject

% this script is used to compute the rotation axis
clear all; close all; clc;

% path information
[MainPath,~] = fileparts(mfilename('fullpath'));
DataPath  = fullfile(MainPath,'Data');

% add functions to the matlab path
addpath(genpath(fullfile(MainPath,'Functions')));

% additional settings
nPP = 81;
BoolGetAxis = false;

%% Get the axis of rotation of the steer for each subject
%ct = 1;
BoolPlot = false;
Folders = {'Classic','EBike'};
% loop over all subjects
for s = 1:nPP
    ppPath = ['pp_' num2str(s)];
    for f = 1:length(Folders)
        % load the data
        OutPathMat = fullfile(DataPath,ppPath,Folders{f});
        filename = fullfile(OutPathMat,'Data_calibrationSteerAxis.mat');
        if exist(filename,'file')
            % load the mat file (processed with ExampleBatch2)
            load(filename,'Callibration');
            timeSpan = [Callibration.t(1) Callibration.t(end)];
            [Rax,n_steer,n_frame] = GetHingeAxis(Callibration.RFrame,Callibration.RSteer,Callibration.t,timeSpan,BoolPlot);
            save(fullfile(OutPathMat,'RotAxis_Steer.mat'),'Rax','n_steer','n_frame');            
        end
    end
end