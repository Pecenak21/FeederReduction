function[struc] = GetCYMData_switchsetting(FilePath)


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
% 1SectionID,2Location,3EqID,4DeviceNumber,5CoordX,6CoordY,7ClosedPhase,8Locked,
% 9RC,10NStatus,11DemandType,12Value1A,13Value2A,14Value1B,15Value2B,16Value1C,17Value2C,
% 18Value1ABC,19Value2ABC,20PhPickup,21GrdPickup,22Alternate,23PhAltPickup,24GrdAltPickup,
% 25FromNodeID,26FaultIndicator,27Automated,28SensorMode,29Strategic,30RestorationMode,
% 31ConnectionStatus,32ByPassOnRestoration,33Reversible

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_SWITCHSETTING.txt'];

if ~exist(FilePath)
    disp('"Switches" structure not populated.');
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
struc.ClosedPhase=y{7}; % Guessing that 'NONE' means that the switch is open
struc.Lock=y{8};
struc.ConnectionStatus=y{31}; % Sounds like this would give the state of the switch ('open' or 'closed') but all are in closed state, which seems unlikely








    