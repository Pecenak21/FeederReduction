function[struc] = GetCYMData_source(FilePath)


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
% [SOURCE]
% FORMAT_SOURCE=1SourceID,2DeviceNumber,3NodeID,4ConnectorIndex,5NetworkID,6DesiredVoltage,7StructureID,8HarmonicEnveloppe,9DemandType,10Value1A,11Value2A,12Value1B,13Value2B,14Value1C,15Value2C,16Value1ABC,17Value2ABC

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_SOURCE.txt'];

if ~exist(FilePath)
    str=[FilePath ' does not exist.'];
    return
else
    fid = fopen(FilePath, 'r');
    y=textscan(fid,strFormat, 'CollectOutput', 0,'delimiter',strDelimiter,'MultipleDelimsAsOne',0);
    fclose(fid);
end
struc.SourceID=y{1};
struc.DeviceNumber=y{2};
struc.NodeID=fnSanitize(y{3});
struc.NetworkID=y{5};





    