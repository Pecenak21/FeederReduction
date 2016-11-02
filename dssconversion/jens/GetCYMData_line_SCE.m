function[struc] = GetCYMData_line(FilePath)


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
% FORMAT_LINE=1ID,2PhaseCondID,3NeutralCondID,4SpacingID,5R1,6R0,7X1,8X0,9B1,10B0,11Amps,12Amps_2,13LockImpedance,14Comments

% HydroOttawa format
% FORMAT_LINE=1ID,2PhaseCondID,3NeutralCondID,4SpacingID,5R1,6R0,7X1,8X0,9B1,10B0,
% 11Amps,12Amps_1,13Amps_2,14Amps_3,15Amps_4,16LockImpedance,17PositiveSequenceShuntConductance,
% 18ZeroSequenceShuntConductance,19Comments

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_LINE.txt'];

if ~exist(FilePath)
    disp('"line" structure not populated');
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
struc.LineID=fnSanitize(y{1});
struc.PhaseCondID=y{2};
struc.NeutralCondID=y{3};
struc.SpacingID=y{4};
struc.R1=y{5};
struc.R0=y{6};
struc.X1=y{7};
struc.X0=y{8};
struc.B1=y{9};
struc.B0=y{10};   