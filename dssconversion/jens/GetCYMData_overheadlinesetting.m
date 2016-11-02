function[struc] = GetCYMData_overheadlinesetting(FilePath)


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
% FORMAT_OVERHEADLINESETTING=1SectionID,2DeviceNumber,3LineCableID,4Length,
% 5ConnectionStatus,6CoordX,7CoordY,8HarmonicModel,9FlowConstraintActive,
% 10FlowConstraintUnit,11MaximumFlow,12SeriesCompensationActive,
% 13MaxReactanceMultiplier,14SeriesCompensationCost

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_OVERHEADLINESETTING.txt'];

if ~exist(FilePath)
    disp('"overheadlinesetting" structure not populated.');
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
struc.ConnectionStatus=y{5};

% struc.CoordX=y{6};
% struc.CoordY=y{7};





    