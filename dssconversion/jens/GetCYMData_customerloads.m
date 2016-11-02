function[struc] = GetCYMData_customerloads(FilePath)


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
% 1SectionID,2DeviceNumber,3LoadType,4CustomerNumber,5CustomerType,
% 6ConnectionStatus,7LockDuringLoadAllocation,8Year,9LoadModelID,
% 10NormalPriority,11EmergencyPriority,12ValueType,13LoadPhase,14Value1,15Value2,
% 16ConnectedKVA,17KWH,18NumberOfCustomer,19CenterTapPercent

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_CUSTOMERLOADS.txt'];

if ~exist(FilePath)
    disp('"customerloads" structure not populated.');
    disp('File not found.');
    FilePath
    struc=[];
    return
else
    fid = fopen(FilePath, 'r');
    y=textscan(fid,strFormat, 'CollectOutput', 0,'delimiter',strDelimiter,'MultipleDelimsAsOne',0);
    fclose(fid);
end
struc.SectionID=fnSanitize(y{1});
% struc.NodeID=y{2};
struc.DeviceNumber=y{2};
struc.LoadType=y{3};
struc.CustomerNumber=y{4};
struc.CustomerType=y{5};
struc.ValueType=y{12};
struc.LoadPhase=y{13};
struc.Value1=y{14};
struc.Value2=y{15};
struc.ConnectedKVA=y{16};
struc.KWH=y{17};
struc.NumberOfCustomer=y{18};










    