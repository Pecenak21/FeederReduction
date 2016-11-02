function[struc struc2] = GetCYMData_LINECONFIGURATION(FilePath)


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
% FORMAT_LINECONFIGURATION=SectionID,DeviceNumber,LineCableID,Length,Overhead,NumberOfCableInParallel,ConnectionStatus,CoordX,CoordY
% 1SectionID,2DeviceNumber,3LineCableID,4Length,5Overhead,6NumberOfCableInParallel,7ConnectionStatus,8CoordX,9CoordY 
% FORMAT_LINECONFIGURATION=SectionID,DeviceNumber,CondID_A,CondID_B,CondID_C,CondID_N,SpacingID,Length,ConnectionStatus,ConductorPosition,CoordX,CoordY
% 1SectionID,2DeviceNumber,3CondID_A,4CondID_B,5CondID_C,6CondID_N,7SpacingID,8Length,9ConnectionStatus,10ConductorPosition,11CoordX,12CoordY

strFormat='%s %s %s %s %s %s %s %s %s';
strFormat2='%s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_LINECONFIGURATION.txt'];

if ~exist(FilePath)
    str=[FilePath ' does not exist.'];
    return
else
    fid = fopen(FilePath, 'r');
    y=textscan(fid,strFormat, 'CollectOutput', 0,'delimiter',strDelimiter,'MultipleDelimsAsOne',0);
    fclose(fid);
end
[flg index]=ismember('FORMAT_LINECONFIGURATION=SectionID',y{1});
if flg==0
    struc.SectionID=fnSanitize(y{1});
    struc.DeviceNumber=y{2};
    struc.LineCableID=y{3};
    struc.Length=y{4};
    struc.Overhead=y{5};
    struc.NumberOfCableInParallel=y{6};
    struc.ConnectionStatus=y{7};
    struc.CoordX=y{8};
    struc.CoordY=y{9};
else
    struc.SectionID=fnSanitize(y{1}(1:index-1));
    struc.DeviceNumber=y{2}(1:index-1);
    struc.LineCableID=y{3}(1:index-1);
    struc.Length=y{4}(1:index-1);
    struc.Overhead=y{5}(1:index-1);
    struc.NumberOfCableInParallel=y{6}(1:index-1);
    struc.ConnectionStatus=y{7}(1:index-1);
    struc.CoordX=y{8}(1:index-1);
    struc.CoordY=y{9}(1:index-1);
    if ~exist(FilePath)
        str=[FilePath ' does not exist.'];
        return
    else
        fid = fopen(FilePath, 'r');
        z=textscan(fid,strFormat2, 'CollectOutput', 0,'delimiter',strDelimiter,'MultipleDelimsAsOne',0,'HeaderLines',index);
        fclose(fid);
        struc2.SectionID=fnSanitize(z{1});
        struc2.DeviceNumber=z{2};
        struc2.CondID_A=z{3};
        struc2.CondID_B=z{4};
        struc2.CondID_C=z{5};
        struc2.CondID_N=z{6};
        struc2.SpacingID=z{7};
        struc2.Length=z{8};
        struc2.ConnectionStatus=z{9};
        struc2.ConductorPosition=z{10};
        struc2.CoordX=z{11};
        struc2.CoordY=z{12};
    end
end







    