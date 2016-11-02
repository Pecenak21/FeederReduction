function[struc] = GetCYMData_shuntcapacitorsetting(FilePath)


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
% 1SectionID,2DeviceNumber,3Location,4Connection,5Phase,6KVAR,7ThreePhaseKVAR,8Losses,
% 9 ThreePhaseLosses,10 KVLN,11 Control,12 ONSetting,13 OFFSetting,14 ShuntCapacitorID,
% 15 ControllingPhase,16 ConnectionStatus,17 SwitchedCapacitorStatus

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_SHUNTCAPACITORSETTING.txt'];

if ~exist(FilePath)
    disp('"Capacitors" structure not populated.');
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
struc.Connection=y{4};
struc.ConnectedPhases=y{5};
struc.kVAR=y{6};
struc.kVLN=y{10}; 






    