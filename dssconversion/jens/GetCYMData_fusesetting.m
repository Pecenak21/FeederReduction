function[struc] = GetCYMData_fusesetting(FilePath)


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
% 1SectionID,2Location,3EqID,4DeviceNumber,5CoordX,6CoordY,7ClosedPhase,8Locked,9RC,
% 10NStatus,11DemandType,12Value1A,13Value2A,14Value1B,15Value2B,16Value1C,17Value2C,18Value1ABC,
% 19Value2ABC,20TCCID,21PhPickup,22GrdPickup,23Alternate,24PhAltPickup,25GrdAltPickup,
% 26FromNodeID,27FaultIndicator,28Strategic,29RestorationMode,30ConnectionStatus,31ByPassOnRestoration,32Reversible

% strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_FUSESETTING.txt'];

if ~exist(FilePath)
    disp('"Fuses" structure not populated');
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
struc.Location=y{2};
struc.EqID=y{3};
struc.DeviceNumber=y{4};
struc.ClosedPhase=y{7}; % Guessing that 'NONE' means that the fuse is open
struc.Lock=y{8};
struc.ConnectionStatus=y{30}; % Sounds like this would give the state of the fuse ('open' or 'closed') but all are in closed state, which seems unlikely








    