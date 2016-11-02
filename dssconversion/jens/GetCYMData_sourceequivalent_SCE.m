function[struc] = GetCYMData_sourceequivalent_SCE(FilePath)


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
% 1NodeID,2Voltage,3PhaseAngle,4PositiveSequenceResistance,5PositiveSequenceReactance,6ZeroSequenceResistance,7ZeroSequenceReactance,8Configuration,9KVLLdesired,10ImpedanceUnit

% Format below for HydroOttawa system:
% 1NodeID,2Voltage,3OperatingAngle1,4OperatingAngle2,5OperatingAngle3,
% 6PositiveSequenceResistance,7PositiveSequenceReactance,8ZeroSequenceResistance,
% 9ZeroSequenceReactance,10Configuration,11OperatingVoltage1,12OperatingVoltage2,
% 13OperatingVoltage3,14ImpedanceUnit

strFormat='%s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_SOURCE EQUIVALENT.txt'];

if ~exist(FilePath)
    disp('"sourceequivalent" structure not populated.');
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
struc.Voltage=y{2};
struc.R1=y{4};
struc.X1=y{5};
struc.R0=y{6};
struc.X0=y{7};
struc.OperatingVoltage1=y{9};
struc.OperatingVoltage2=y{9};
struc.OperatingVoltage3=y{9};
struc.ImpedanceUnit=y{10};






    