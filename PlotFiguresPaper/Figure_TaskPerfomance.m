%% Figure 3

% path information
[TempPath,~,~] = fileparts(mfilename('fullpath'));
MainPath = TempPath(1:end-17);
FileTaskInfo = fullfile(MainPath,'Data','TaskInformation.xlsx');

% add functions to the matlab path
addpath(genpath(fullfile(MainPath,'Functions')));

% read the excel file with quantitative information
[ShoulderCheckInfo] = GetShoulderCheckInfo(FileTaskInfo);

% open a figure
figure();
set(gcf,'Position',[634   360   794   360]);

% plot jongeren (no errors)
subplot(1,2,1)
pie([10 0 0 0 0]);

% plot ouderen
subplot(1,2,2)
Zadel = sum(ShoulderCheckInfo.DatOlder(:,1)==1 & ShoulderCheckInfo.DatOlder(:,3)~=1);
Kegel = sum(ShoulderCheckInfo.DatOlder(:,2)==1 & ShoulderCheckInfo.DatOlder(:,3)~=1);
Voet = sum(ShoulderCheckInfo.DatOlder(:,3)==1);
BuitLijn = sum(ShoulderCheckInfo.DatOlder(:,4)==1 & ShoulderCheckInfo.DatOlder(:,3)~=1);
AllesOK = sum(sum(ShoulderCheckInfo.DatOlder(:,1:4),2) == 0);
pie([AllesOK BuitLijn Voet Kegel Zadel]);


