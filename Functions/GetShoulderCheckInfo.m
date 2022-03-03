function [ShoulderCheckInfo] = GetShoulderCheckInfo(filename)
%GetShoulderCheckInfo Reads all the information related to the shoulder
%check task from the excel table (that was logged by the test leaders)

opts = spreadsheetImportOptions('NumVariables',16);
opts.DataRange = 'K4:Z44';
opts.Sheet = 'Info-Oud';
DatOlder = readcell(filename,opts);
DatOlder = cell2mat(DatOlder );

opts = spreadsheetImportOptions('NumVariables',16);
opts.DataRange = 'K3:Z3';
opts.Sheet = 'Info-Oud';
HeadersOlder = readcell(filename,opts);

opts = spreadsheetImportOptions('NumVariables',16);
opts.DataRange = 'M5:AB44';
opts.Sheet = 'Info-Jong';
DatYoung= readcell(filename,opts);
DatYoung = cell2mat(DatYoung );

opts = spreadsheetImportOptions('NumVariables',16);
opts.DataRange = 'M4:AB4';
opts.Sheet = 'Info-Jong';
HeadersYoung = readcell(filename,opts);

opts = spreadsheetImportOptions('NumVariables',1);
opts.DataRange = 'A5:A44';
opts.Sheet = 'Info-Jong';
ppIDYoung= readcell(filename,opts);
ppIDYoung = cell2mat(ppIDYoung);

opts = spreadsheetImportOptions('NumVariables',1);
opts.DataRange = 'A4:A44';
opts.Sheet = 'Info-Oud';
ppIDOld= readcell(filename,opts);
ppIDOld = cell2mat(ppIDOld);


ShoulderCheckInfo.DatYoung = DatYoung;
ShoulderCheckInfo.DatOlder = DatOlder;
ShoulderCheckInfo.HeadersOlder = HeadersOlder;
ShoulderCheckInfo.HeadersYoung = HeadersYoung;
ShoulderCheckInfo.ppIDYoung = ppIDYoung;
ShoulderCheckInfo.ppIDOld = ppIDOld;
end

