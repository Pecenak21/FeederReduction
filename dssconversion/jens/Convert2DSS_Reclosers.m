function rc = Convert2DSS_Reclosers(rcs)
% rc = Convert2DSS_Reclosers(rcs)
%
% PURPOSE : Converts recloser from object to OpenDSS
%
%
% INPUT :   rcs: struct that has all entries (e.g., o.R1)
%
% OUTPUT :  OpenDSS recloser object
%           Not documented in OpenDSS manual V7.6(November 2012)
%           Information in help file of OpenDSS V7.6.2.1:

rc(length(rcs)) = dssrecloser;


rc.Name = rcs.DeviceNumber;
rc.MonitoredObj = ['Line.' rcs.SectionID];

if isfield(rcs,'FastPhaseCount') && isfield(rcs,'SlowPhaseCount') ...
        && isfield(rcs,'SlowGroundTimeAddSec') && isfield(rcs,'FastGroundTimeAddSec') ... 
        && isfield(rcs,'SlowPhaseTimeAddSec') && isfield(rcs,'FastPhaseTimeAddSec')
    % All this info is in SynerGEE, but does not seem to be in CYME
    if(~rcs.NearFromNode)
        rc.SwitchedTerm = 2;
    end
    rc.NumFast = rcs.FastPhaseCount;
    rc.Shots = rcs.FastPhaseCount + rcs.SlowPhaseCount;
    rc.TDGrDelayed = rcs.SlowGroundTimeAddSec;
    rc.TDGrFast = rcs.FastGroundTimeAddSec;
    rc.TDPhDelayed = rcs.SlowPhaseTimeAddSec;
    rc.TDPhFast = rcs.FastPhaseTimeAddSec;
else
    % Model as a switch if no detailed recloser info is available
    if strcmp(rcs.ClosedPhase,'NONE')
        rc.Action = 'Open';
    else
        rc.Action = 'Close';
    end    
    
end

