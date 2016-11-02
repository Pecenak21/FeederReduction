function[struc] = GetCYMData_cable(FilePath)


% GetPSSEData_Bus(FilePath)
%
% PURPOSE : get content of bus parameters in file
%
%
% INPUT :   FilePath: directory and file name
%           strFormat: format of data in file, 
%               e.g. '%s %f %f %f' 1st column is string , 2nd through 4th
%               columns are floating point numbers
%           strDelimiter: delimiter that separates columns in data file
%           flgHeader: 0 -> no header, 1 -> header
%           strDelimiter_Header: delimiter that separates header columns 
%
% OUTPUT :  y_First: First column in cell array
%           y: Complete cell array
%           y_Header: Header array (if flgHeader is set)
%           struc: a structure with header names and respective column content under each header name (if flgHeader is set)
%
% FORMAT_CONCENTRICNEUTRALCABLE=1ID,2R1,3R0,4X1,5X0,6B1,7B0,8Amps,9Amps_1,
% 10Amps_2,11Amps_3,12Amps_4,13WithstandRating,14NumberOfPhases,15FailRate,
% 16TmpFailRate,17MajorRepairTime,18MinorRepairTime,19MajorFailureProportion,
% 20PhaseCondID,21NeutralCondID,22NeutralWires,23InsulationType,
% 24InsulationDielectricConst,25InsulationDiameter,26Dab,27Dbc,28Dac,29RatedLevel,
% 30ContinuousTemperatureRating,31SCCurrentTemperatureRating,
% 32PositiveSequenceShuntConductance,33ZeroSequenceShuntConductance,
% 34LockImpedance,35Comments

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_CONCENTRICNEUTRALCABLE.txt'];

if ~exist(FilePath)
    disp('"cable" structure not populated');
    disp('File not found.');
    FilePath
    struc=[];
    return
else
    fid = fopen(FilePath, 'r');
    y=textscan(fid,strFormat, 'CollectOutput', 0,'delimiter',strDelimiter,'MultipleDelimsAsOne',0);
    fclose(fid);
end
%y{1}
% struc.Info=y;
struc.CableID=y{1};
%y{2}
struc.R1=y{2};
%y{3}
struc.R0=y{3};
struc.X1=y{4};
struc.X0=y{5};
struc.B1=y{6};
struc.B0=y{7};






    