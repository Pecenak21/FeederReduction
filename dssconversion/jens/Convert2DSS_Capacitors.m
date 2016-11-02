function [cp capcon] = Convert2DSS_Capacitors(cps)
% cp = Convert2DSS_Capacitors(cps)
%
% PURPOSE : Converts capacitor from object to OpenDSS, creates capacitor
%           objects and capacitor control objects
%
%
% INPUT :   cps: struct that has all entries (e.g., o.R1)
%
% OUTPUT :  OpenDSS Capacitor object
%           see page 110 in OpenDSS manual V7.6 (November 2012)
%           OpenDSS CAPControl object
%           see page 143 in OpenDSS manual V7.6 (November 2012)

% Capacitor
cp(length(cps)) = dsscapacitor;
% corresponding capacitor control object
capcon(length(cps)) = dsscapcontrol;

for i = 1:length(cps)
    % Create capacitor object
    cp(i).Name = cps(i).DeviceNumber; % currently not used, .UniqueDeviceId in SynerGEE
%     buses = regexp(cps(i).SectionID, '_', 'split');
%     cp(i).Bus1 = cps(i).FromNodeID;
    cp(i).Bus1 = cps(i).ToNodeID; % capacitor shunt connected to 'to' bus
    cp(i).Phases = cps(i).Phase; % .ConnectedPhases in SynerGEE
    cp(i).Kv = str2double(cps(i).KVLN); % .RatedKv in SynerGEE
    ThreePhaseKVAR=str2double(cps(i).ThreePhaseKVAR);
    if ThreePhaseKVAR==0 % sometimes 'ThreePhaseKVAR' is erroneously specified as '0'
        KVAR=str2double(cps(i).KVAR);
        ThreePhaseKVAR=3*KVAR;
    end
    cp(i).Kvar = ThreePhaseKVAR; % this is the kVAR rating for all phases combined
    cp(i).Conn = connClean(cps(i).Connection);
    if isfield(cps,'CapacitorIsOn')
        cp(i).states = cps(i).CapacitorIsOn; % Not sure where this is specified in CYME. Leaving OpenDSS default, which is 'on'
    else
        cp(i).states = 1; % default, specifying it makes this parameter show up in the write file
    end
    % Create CAPControl object
    capcon(i).Name = ['capctrl_' cps(i).DeviceNumber];
	capcon(i).Element = ['line.' cps(i).SectionID];
	capcon(i).Capacitor = cps(i).DeviceNumber;
    if isfield(cps,'PrimaryControlMode') % SynerGEE parameter, not sure how to get this from CYME
        capcon(i).Type = cps(i).PrimaryControlMode;
        switch capcon(i).Type
            case 'voltage'
                capcon(i).PTPhase = cps(i).MeteringPhase;
                capcon(i).OFFsetting = nonzeros([cps(i).Module1CapSwitchTripValue cps(i).Module2CapSwitchTripValue cps(i).Module3CapSwitchTripValue]);
                capcon(i).ONsetting = nonzeros([cps(i).Module1CapSwitchCloseValue cps(i).Module2CapSwitchCloseValue cps(i).Module3CapSwitchCloseValue]);
            case 'current'
                capcon(i).CTPhase = cps(i).MeteringPhase;
                capcon(i).OFFsetting = nonzeros([cps(i).Module1CapSwitchTripValue cps(i).Module2CapSwitchTripValue cps(i).Module3CapSwitchTripValue]);
                capcon(i).ONsetting = nonzeros([cps(i).Module1CapSwitchCloseValue cps(i).Module2CapSwitchCloseValue cps(i).Module3CapSwitchCloseValue]);
            case 'kvar' % Not tested for 'kvar'
                capcon(i).OFFsetting = min(nonzeros([cps(i).Module1CapSwitchTripValue cps(i).Module2CapSwitchTripValue cps(i).Module3CapSwitchTripValue]));
                capcon(i).ONsetting = max(nonzeros([cps(i).Module1CapSwitchCloseValue cps(i).Module2CapSwitchCloseValue cps(i).Module3CapSwitchCloseValue]));
        end
        capcon(i).VoltOverride = cps(i).CapVoltageOverrideActive;
        capcon(i).Vmin = cps(i).CapVoltageOverrideSetting - cps(i).CapVoltageOverrideBandwidth/2;
        capcon(i).Vmax = cps(i).CapVoltageOverrideSetting + cps(i).CapVoltageOverrideBandwidth/2;
        capcon(i).PTRatio = cps(i).CapacitorPTRatio;
        capcon(i).CTRatio = cps(i).CapacitorCTRating;
    end


end
    
% obj.defaults = struct('Name',{''}, ...
% 		'Like', {''}, ...
% 		'Phases',3, ...
% 		'Bus1','', ...
% 		'Bus2','', ...
% 		'Numsteps',1, ... numsteps normally comes after kvar in dss notation, however based on how we write output files, we wouldn't be able to specify different amounts for different steps that way; instead we change the order, and lose the splitting behavior that opendss offers.
% 		'Kvar',1200, ... the documentation says the default is 600, but when you actually make one, you get 1200 by default in the version I have here.
% 		'Kv',12.47, ...
% 		'Conn','wye', ...
% 		'Cmatrix',{[]}, ...
% 		'Cuf',20.47, ... Per phase; default based on default kvar, kv, and basefreq values
% 		'R',0, ...
% 		'XL',0, ...
% 		'Harm',0, ...
% 		'states',1, ...
% 		'Normamps',{[]}, ...
% 		'Emergamps',{[]}, ...
% 		'Faultrate',0.0005, ...
% 		'Pctperm',100, ...
% 		'Basefreq',60);

% 	obj.defaults = struct(...
% 		'Name',{''}, ...
% 		'Like','',...
% 		'Element','', ... %full name (e.g. 'Line.line1')
% 		'Capacitor','', ... %not full name (e.g. 'cap1')
% 		'Type','current', ... %one of {'current','voltage','kvar','pf','time'}  % I had to read the code to find the default so be glad you have it
% 		'CTPhase',1, ... % phase number; for delta/ll connection use the first of the two phases (so 1 for 1-2, 2 for 2-3, 3 for 3-1)
% 		'CTRatio',60, ... % ratio between line amps and control amps
% 		'DeadTime',300, ... % time after OFF before ON allowed (seconds)
% 		'Delay',15, ... %time between control being armed and turning on (during which time, action can be cancelled) in seconds
% 		'DelayOFF',15, ... % same but for turning off
% 		'EventLog','YES', ... % YES/true or NO/false on whether to log actions
% 		'OFFsetting',200, ... % value of control variable at which to switch off
% 		'ONsetting',300, ... % value of control variable at which to switch on
% 		'PTPhase',1, ... %phase to monitor for voltage control
% 		'PTRatio',60, ... % ratio between monitored voltage and control voltage
% 		'terminal',1, ... % which terminal of the element at which to monitor
% 		'VBus','', ... % bus to use for voltage override (below).  Defaults to bus of terminal.
% 		'Vmax',126, ...% max voltage (after PT ratio) regardless of other controls
% 		'Vmin',115, ...% min voltage (after PT ratio) regardless of other controls
% 		'VoltOverride','NO'); % whether to use Vmax/Vmin
% 	obj.fieldnames = fieldnames(obj.defaults);
% 	% set the data structure to be a blank version of the defaults
% 	% structure
% 	obj.data = obj.defaults;
% 	obj.data(1) = []; obj.data(1).Name = ''; obj.data(1).Like = '';

% If you wish a series capacitor, simply specify a second bus connection.
% If you wish an ungrounded wye capacitor, set all the second terminal conductors to an empty
% node on the first terminal bus, e.g.:
% Bus1=B1 bus2 = B1.4.4.4 ! for a 3-phase capacitor
% Of course, any other connection is possibly by explicitly specifying the nodes.
% While the Capacitor object may be treated as if it were a single capacitor bank, it is actually
% implemented as a multistep tuned filter bank. Many of the properties can be specified as arrays.
% See Numsteps property.
% If you wish to represent a filter bank, specify either the XL or Harm property. When a Capcontrol
% object switches a capacitor, it does so by incrementing or decrementing the active step.
% Properties, in order, are:
% Bus1= Definition for the connection of the first bus. When this is set, Bus2 is set to the
% same bus name, except with all terminals connected to node 0 (ground reference).
% Set Bus 2 at some time later if you wish a different connection.
% Bus2= Bus connection for second terminal. Must always be specified after Bus1 for series
% capacitor. Not necessary to specify for delta or grd?wye shunt capacitor. Must be
% specified to achieve an ungrounded neutral point connection.
% Phases= Number of phases. Default is 3.
% Kvar= Most common of three ways to define a power capacitor. Rated kvar at rated kV,
% total of all phases. Each phase is assumed equal. Normamps and Emergamps
% automatically computed. Default is 600.0 kvar. Total kvar, if one step, or ARRAY of
% kvar ratings for each step. Evenly divided among phases. See rules for NUMSTEPS.
% Kv= Rated kV of the capacitor (not necessarily same as bus rating). For Phases=2 or
% Phases=3, enter line?to?line (phase?to?phase) rated voltage. For all other numbers
% of phases, enter actual can rating. (For Delta connection this is always line?to?line
% rated voltage). Default is 12.47 kV.
% Conn= Connection of bank. One of {wye | ln} for wye connected banks or {delta | ll} for
% delta (line?line) connected banks. Default is wye (or straight?through for series
% capacitor).
% Cmatrix= Alternate method of defining a capacitor bank. Enter nodal capacitance matrix in
% ?f. Can be used to define either series or shunt banks. Form should be:
% Cmatrix=[c11 | ?c21 c22 | ?c31 ?c32 c33]
% All steps are assumed the same if this property is used.
% Cuf= Alternate method of defining a capacitor bank. Enter a value or ARRAY of values for
% C in ?f. ARRAY of Capacitance, each phase, for each step, microfarads. See Rules for
% NumSteps.
% R = ARRAY of series resistance in each phase (line), ohms. Default is 0.0
% XL = ARRAY of series inductive reactance(s) in each phase (line) for filter, ohms at base
% frequency. Use this OR "h" property to define filter. Default is 0.0.
% Harm = ARRAY of harmonics to which each step is tuned. Zero is interpreted as meaning
% zero reactance (no filter). Default is zero.
% Numsteps = Number of steps in this capacitor bank. Default = 1. Forces reallocation of the
% capacitance, reactor, and states array. Rules: If this property was previously =1, the
% value in the kvar property is divided equally among the steps. The kvar property
% does not need to be reset if that is correct. If the Cuf or Cmatrix property was used
% previously, all steps are set to the value of the first step. The states property is set
% to all steps on. All filter steps are initially set to the same harmonic. If this property
% was previously >1, the arrays are reallocated, but no values are altered. You must
% SUBSEQUENTLY assign all array properties.
% states = ARRAY of integers {1|0} states representing the state of each step (on|off). Defaults
% to 1 when reallocated (on). Capcontrol will modify this array as it turns steps on or
% off.
% Normamps= Normal current rating. Automatically computed if kvar is specified. Otherwise, you
% need to specify if you wish to use it.
% Emergamps= Overload rating. Defaults to 135% of Normamps.
% Faultrate= Annual failure rate. Failure events per year. Default is 0.0005.
% Pctperm= Percent of faults that are permanent. Default is 100.0.
% Basefreq= Base frequency, Hz. Default is 60.0
% Like= Name of another Capacitor object on which to base this one.
