function[struc] = GetCYMData_loads(FilePath,flgTest)


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
% Format: 1SectionID,2DeviceNumber,3LoadType,4Connection,5Location

if nargin<2
    flgTest=0;
end

strFormat='%s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);



if flgTest
    FilePath=[FileName '_LOADS_small.txt'];
else
    FilePath=[FileName '_LOADS.txt'];
end

if ~exist(FilePath)
    disp('"Loads" structure not populated.');
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
struc.LoadType=y{3};





    