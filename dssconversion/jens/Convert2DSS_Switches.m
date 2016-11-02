function sw = Convert2DSS_Switches(sws,useuid)
% sw = Convert2DSS_Switches(sws,useuid)
%
% PURPOSE : Converts switch from object to OpenDSS
%
%
% INPUT :   sws: struct that has all entries (e.g., o.R1)
%           useuid: switch has unique device ID
%
% OUTPUT :  OpenDSS switch object
%           Not documented in OpenDSS manual V7.6(November 2012)
%           Information in help file of OpenDSS V7.6.2.1:

sw(length(sws)) = dssswtcontrol;


% Sometimes the "unique device ID" field isn't actually unique.
% When it is, we use it, otherwise we use the section ID:
if(useuid)
    sw.Name = sws.DeviceNumber; % sws.UniqueDeviceId;
else
    sw.Name = sws.SectionID;
end

% object of the switch in this case is the line/section
sw.SwitchedObj = ['Line.' sws.SectionID];

if strcmp(sws.ClosedPhase,'NONE')
    sw.Action = 'Open';
else
    sw.Action = 'Close';
end

%TODO: check what kind of data type Lock field receives
% sw.Lock = logical(sws.SwitchIsTie);

% if(~sws.NearFromNode)
%     sw.SwitchedTerm = 2;
% end

