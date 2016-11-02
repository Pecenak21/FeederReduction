function[struc] = GetCYMData_undergroundlinesetting(FilePath)


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
% FORMAT_UNDERGROUNDLINESETTING=1SectionID,2DeviceNumber,3LineCableID,
% 4Length,5NumberOfCableInParallel,6Amps,7Amps_1,8Amps_2,9Amps_3,10Amps_4,
% 11ConnectionStatus,12CoordX,13CoordY,14HarmonicModel,15FlowConstraintActive,
% 16FlowConstraintUnit,17MaximumFlow 

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_UNDERGROUNDLINESETTING.txt'];

if ~exist(FilePath)
    disp('"undergroundlinesetting" structure not populated.');
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
struc.DeviceNumber=y{2};
struc.LineCableID=y{3};
struc.Length=y{4};
struc.NumberOfCableInParallel=y{5};
struc.ConnectionStatus=y{11};
% struc.CoordX=y{12};
% struc.CoordY=y{13};





    