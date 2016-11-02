function[struc] = GetCYMData_node(FilePath)


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
% Format:
% 1NodeID,2CoordX,3CoordY,4TagText,5TagProperties,6TagDeltaX,7TagDeltaY,8TagAngle,
% 9TagAlignment,10TagBorder,11TagBackground,12TagTextColor,13TagBorderColor,
% 14TagBackgroundColor,15TagLocation,16TagFont,17TagTextSize,18TagOffset,19ZoneID,
% 20ConnectedEquipType,21ExposedCircuitType,22BusGap,23WorkingDistance,
% 24UseUserDefinedFaultCurrent,25UserDefinedFaultCurrent,26OpeningTimeMode,
% 27UserDefinedOpeningTime,28OverrideLFVoltageLimit,29HighVoltageLimit,
% 30LowVoltageLimit,31LoadSheddingActive,32MaximumLoadShed,33ShedLoadCost,34UserDefinedBaseVoltage

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_NODE.txt'];

if ~exist(FilePath)
    disp('"Node" structure not populated');
    disp('File not found.');
    FilePath
    struc=[];
    return
else
    fid = fopen(FilePath, 'r');
    y=textscan(fid,strFormat, 'CollectOutput', 0,'delimiter',strDelimiter,'MultipleDelimsAsOne',0);
    fclose(fid);
end
struc.NodeID=fnSanitize(y{1});
struc.CoordX=y{2};
struc.CoordY=y{3};







    