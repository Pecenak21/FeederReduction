function output = dssconversion(circuitdata, linecodedata, circuitname)
% This conversion engine works for data from SynerGee data (in Excel format) only.
% Modification is needed if other source of data is used.
% Inputs:
%			circuitdata: circuit's data or path to excel file
%			linecodedata: linecode data or path to linecode excel file

%			circuitname: (optional, ignored). Default: 'newcircuit'

% Process inputs
if isstruct(circuitdata)
	d = circuitdata;
elseif isa(circuitdata, 'char')
	d = excel2obj(circuitdata);
else
	error('Invalid circuitdata type');
end

if isstruct(linecodedata)
	glc = linecodedata;
elseif isa(linecodedata, 'char')
	fprintf('Reading linecodes... '); t_ = tic;
	glc = excel2obj(linecodedata);
	glc = glc.LineCode;
	toc(t_);
else
	error('Invalid linecodedata type');
end

% keep track of fieldnames we handle
% when we're all done we'll report unhandled fields to the user
fns_handled = {};

%% Linecodes


%% added by changfu to get rid of bracket in the new linecode got from SDG&E
for i = 1:length(glc);
    glc(i).ConductorId = strrep(glc(i).ConductorId,'(','');
    glc(i).ConductorId = strrep(glc(i).ConductorId,')','');
end

fprintf('Converting linecodes... '); t_ = tic;
lc(length(glc)) = dsslinecode;
for i = 1:length(glc)
    
	lc(i).Name = glc(i).ConductorId;
    
%     i
%     glc(i).ConductorId
%     lc(i).Name

	lc(i).R1 = glc(i).R1;
	lc(i).X1 = glc(i).X1;
	lc(i).R0 = glc(i).R0;
	lc(i).X0 = glc(i).X0;
	lc(i).Units = 'kft';
end
toc(t_);

%% Sections / lines
fprintf('Converting sections... '); t_ = tic;
s = d.InstSection;
l(length(s)) = dssline;
% make phase mapping easier:
phasemap.(d.SAI_Control.Phase1) = '.1';
phasemap.(d.SAI_Control.Phase2) = '.2';
phasemap.(d.SAI_Control.Phase3) = '.3';
phasemap.N = '.0';
if(any(strcmpi(d.SAI_Control.LengthUnits,{'English','English2'})))
	lineunits = 'Ft';
elseif(strcmpi(d.SAI_Control.LengthUnits,'Metric'))
	lineunits = 'meter';
end
for i = 1:length( s )
	l(i).Name = s(i).SectionId;
    % setting the linecode
    % 1. search for the linecode id
    
    % tmp{1} = dataclean(s(i).PhaseConductorId,'name');
    
    %%
    
    % added by changfu to deal with bracket in the new linecode from SDG&E
    % refer to line 37 in this function
    s(i).PhaseConductorId = strrep(s(i).PhaseConductorId,'(','');
    s(i).PhaseConductorId = strrep(s(i).PhaseConductorId,')','');
 
    
    [x idx] = ismember(dataclean(s(i).PhaseConductorId,'name'), get(lc(:),'Name'));
    if(~all(x))
        error('dssConversion:missingLinecode', 'can''t find line code named %s', dataclean(s(i).PhaseConductorId,'name'));
    end
    % 2. set linecode of the line to be that linecode object
	l(i).LineCode = lc(idx);

	%only report a neutral phase if we're going to ground it
	secphases = s(i).SectionPhases;
	secphases(secphases==' ') = [];
	if(~s(i).NeutIsGrounded)
		secphases(secphases=='N') = [];
	end
	l(i).Phases = length(secphases);
	% convert units to the units we'll be using for everything else (actually just length right now)
    l(i).Units = lineunits;
	% do length
	l(i).Length = s(i).SectionLength_MUL; %MUL: medium unit length (feet for english or english2 units, and meters for metric units)
	secphstr = '';
	for phase_idx = 1:length(secphases)
		secphstr = [secphstr phasemap.(secphases(phase_idx))];
	end
	l(i).bus1 = [s(i).FromNodeId secphstr];
	l(i).bus2 = [s(i).ToNodeId secphstr];
	if(~s(i).IdenticalPhaseConductors)
		warning('dssconversion:unimplemented','Use of different phases for different conductors is not implemented');
    end
end
fns_handled{end+1} = 'InstSection';
toc(t_);

%% Loads
fprintf('Converting loads... '); t_ = tic;
% get loads in nicer struct format
ldi = d.Loads;
if(length(ldi)==1 && iscell(ldi.SectionId)), ldi = structconv(ldi); end
ld = dssload(); ld(length(ldi)).Name = '';
idx = 1;
phases = [[ldi.Phase1Kw]' [ldi.Phase2Kw]' [ldi.Phase3Kw]'];
[tmpvar sectionind] = ismember({ldi.SectionId},{d.InstSection.SectionId});

for i = 1:length( ldi )

	% get model parameters from the line
	ldsec = d.InstSection(sectionind(i));
	if(ldi(i).IsSpotLoad)
		zipv = [ldsec.PercentSpotLoadConstImpedance, ldsec.PercentSpotLoadConstCurrent];
	else
		zipv = [ldsec.PercentDistLoadConstImpedance, ldsec.PercentDistLoadConstCurrent];
	end
	zipv = [zipv 100-sum(zipv)]/100; % add the power part and convert from percent to fraction
	zipv = [zipv zipv .5]; % assume real and reactive parts same, set cutoff to 0.5pu (may be way too low?
	% for now we ignore the Kwh field for each phase, since it's not clear
	% what that would give to the model, and they're all zero anyway.
	
	% we're currently opting to model all loads as spot loads at the 'to'
	% node unless expicitly specified as 'F' (from node), which doesn't
	% actually happen in any of our models.
	if(ldi(i).LocationToModelSpotLoads == 'F')
		busid = ldsec.FromNodeId;
	else
		busid = ldsec.ToNodeId;
	end
	if(ldi(i).IsSpotLoad || true)
		ld(idx).bus1 = busid;
		% get a list of the power on each of the phases
		kw = phases(i,phases(i,:)>0);
		% Right now I'm assuming only LN connected loads.  Some loads might
		% be LL, in which case that would be weird, but our current model
		% specifies all the sections as having a neutral connection, so we
		% don't need to worry about that here.
		if(length(kw)<3 || length(unique(kw))>1) %if we have unequal power on the phases or less than three phases, we model power use with one load per phase
			for ph=find(phases(i,:)>0)
				ld(idx).Name = ['l' busid '_' d.SAI_Control.(sprintf('Phase%i',ph))];
				ld(idx).bus1 = sprintf('%s.%i',busid,ph);
				ld(idx).Phases = 1;
				ld(idx).NumCust = ldi(i).(sprintf('Phase%iCustomers',ph));
				ld(idx).Kw = phases(i,ph);
				ld(idx).Kvar = ldi(i).(sprintf('Phase%iKvar',ph));
				ld(idx).Model = 8; %ZIP load
				ld(idx).ZIPV = zipv;
				idx = idx+1;
			end
		else % if we have three phases with equal power, we can use a single 3-phase load
			ld(idx).Name = ['l' busid];
			ld(idx).Phases = 3;
			ld(idx).NumCust = ldi(i).Phase1Customers + ldi(i).Phase2Customers + ldi(i).Phase3Customers;
			ld(idx).Kw = sum(kw);
			ld(idx).Kvar = ldi(i).Phase1Kvar+ldi(i).Phase2Kvar+ldi(i).Phase3Kvar;
			ld(idx).Model = 8; %ZIP load
			ld(idx).ZIPV = zipv;
			idx = idx + 1;
		end
	else
		% The treatment above is actually correct only for spot loads.
		% Non spot loads are indicated to be modeled at 'center' which is
		% not what we're doing right now.
	end
end
fns_handled{end+1} = 'Loads';
toc(t_);

%% Capacitors and CapControls
% for each capacitor, a capacitor controller is created based on given
% information.
if(isfield(d,'InstCapacitors'))
	fprintf('Converting capacitors... '); t_ = tic;
	cps = d.InstCapacitors;
    
    [tmpvar sectionind] = ismember({cps.SectionId},{d.InstSection.SectionId});

% 
% for i = 1:length( ldi )
% 
% 	% get model parameters from the line
% 	ldsec = d.InstSection(sectionind(i));


	% capacitor
	cp(length(cps)) = dsscapacitor;

	% corresponding capacitor control object
	capcon(length(cps)) = dsscapcontrol;

	for i = 1:length(cps)
        
        Idsec = d.InstSection(sectionind(i));
		% set up capacitor
		cp(i).Name = cps(i).UniqueDeviceId;
		% buses = regexp(cps(i).SectionId, '_', 'split');
        
		% cp(i).Bus1 = buses{2};
        cp(i).Bus1 = Idsec.ToNodeId;
		cp(i).Phases = cps(i).ConnectedPhases;
		cp(i).Kv = cps(i).RatedKv;
		
		% get fixed kvar rating
		fixedKvar = sum([cps(i).FixedKvarPhase1 cps(i).FixedKvarPhase2 cps(i).FixedKvarPhase3]);
		% get switchable kvar ratings
		kvar = cp(i).Phases*[cps(i).Module1KvarPerPhase cps(i).Module2KvarPerPhase cps(i).Module3KvarPerPhase];
		ind_z = find(kvar==0,1);
		kvar(ind_z:end) = [];
		% set one or the other of them as appropriate
		if fixedKvar > 0
			if(~isempty(kvar)), warning('dssconversion:capacitor:multiplekvar','Capacitor kvar was specified both as fixed and as switchable, which is invalid'); end
			cp(i).Kvar = fixedKvar;
        elseif ~isempty(kvar)
			if(length(kvar)>1), cp(i).Numsteps=length(kvar); end
			cp(i).Kvar = kvar;
		else
			error('DSSConversion:Capcitor:kvar','Invalid kvar input');
		end
		
		cp(i).Conn = cps(i).ConnectionType;
		cp(i).states = cps(i).CapacitorIsOn;

		% set up capControl
		capcon(i).Name = ['capctrl_' cps(i).UniqueDeviceId];
		capcon(i).Element = ['line.' cps(i).SectionId];
		capcon(i).Capacitor = cps(i).UniqueDeviceId;
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
	toc(t_);
	fns_handled{end+1} = 'InstCapacitors';
end


%% Large customers as loads
if(isfield(d,'InstLargeCust'))
	fprintf('Converting large customer loads... '); t_ = tic;
	ldi = d.InstLargeCust;
	if(length(ldi)==1 && iscell(ldi.SectionId)), ldi = structconv(ldi); end
	ld(idx+length(ldi)-1).Name = ''; % allocate another chunk of data
	
	% initialize generators
	genId = find([d.InstLargeCust.GenTotalKw]>0);
	if sum(genId) > 0
		ge(sum(genId)) = dssgenerator;
		
		% intialize generator classes
		genClasses = {'unknown'};
		%new class
		nc = unique(lower({d.InstLargeCust(genId).GenCustClass}));
		genClasses = [genClasses setdiff(nc,genClasses)];
	end
	
	for i = 1:length( ldi )
		ld(idx).Name = ldi(i).UniqueDeviceId;

		% look up ZIP values for the load
		zipv = [ldi(i).LoadPctConstZ, ldi(i).LoadPctConstI];
		zipv = [zipv ldi(i).LoadPct-sum(zipv)]/ldi(i).LoadPct; % add the power part and convert from percent to fraction
		zipv = [zipv zipv .5]; % assume real and reactive parts same, set cutoff to 0.5pu (may be way too low?

		% look up the 'to' or 'from' node of the section
		busid = strcmp(ldi(i).SectionId,{d.InstSection.SectionId});
		ldsec = d.InstSection(busid);

		if(ldi(i).Location == 'F') % indicates "From"
			busid = ldsec.FromNodeId;
		else % would be 'T' for "to" or 'c' for "center", so we use toNode
			busid = ldsec.ToNodeId;
		end
		ld(idx).bus1 = busid;
		kw = [ldi(i).LoadPhase1Kw ldi(i).LoadPhase2Kw ldi(i).LoadPhase2Kw];
		if(sum(kw>0)<3 || length(unique(kw))>1), warning('dssconversion:Unimplemented','Unbalanced large customers are not currently implemented'); end
		ld(idx).Phases = sum( kw>0 );
		ld(idx).NumCust = sum( [ldi(i).Phase1Customers ldi(i).Phase2Customers ldi(i).Phase3Customers]);
		ld(idx).Kw = sum( kw );
		ld(idx).Kvar = sum( [ldi(i).LoadPhase1Kvar ldi(i).LoadPhase2Kvar ldi(i).LoadPhase2Kvar]);
		ld(idx).Model = 8; %ZIP load
		ld(idx).ZIPV = zipv;

		% add generator if customer is generating electricity
		if(ldi(i).GenTotalKw>0)
			if(~exist('gid','var'))
				gid = 1;
			else
				gid = gid + 1;
			end
			ge(gid).Name = ldi(i).UniqueDeviceId; 
			ge(gid).Kw = ldi(i).GenTotalKw * ldi(i).GenPct/100;
			ge(gid).Kvar = ldi(i).GenTotalKvar;
			ge(gid).bus1 = busid;
			
			% PV generator
			if ~isempty(strfind('pv',lower(ldi(i).GenCustClass)))
				ge(gid).Model = 3;
			end

			% set class. Class 0 is 'unknown' (our convention).
			[val c] = ismember(lower(ldi(i).GenCustClass),genClasses);
			ge(gid).Class = c - 1;
			
			if strcmp(ldi(i).GenStatus,'P') %powered
				ge(gid).Enabled = 'yes';
			elseif strcmp(ldi(i).GenStatus,'O') %offline
				ge(gid).Enabled = 'no';
			else
				warning('dssconversion:generatorSetup','check generator status again. You might need to change the codes according to your data on how they specify generator status. Set to YES by default.');
				ge(gid).Enabled = 'yes';
			end
		end

		idx = idx + 1;
	end
	toc(t_);
	fns_handled{end+1} = 'InstLargeCust';
end
%% Switches
if(isfield(d,'InstSwitches'))
	fprintf('Converting switches... '); t_ = tic;
	sws = d.InstSwitches;
	sw(length(sws)) = dssswtcontrol;

	if(any([sws.IsAutomaticSwitch]))
		warning('dssconversion:Unimplemented','Automated switching is currently not implemented; if this feature is important to your project, please add that to the code!');
	end
	if(any([sws.SwitchIsTie]))
		warning('dssconversion:SwitchIsTie:NeedToCheck','Implemented SwitchIsTie as "lock" parameter of switch (best guess). You may need to double check and see if this is what you want in your circuit!');
	end
	% Determine which unique ids are actually unique (used for naming)
	% first get them all
	useuid = {sws.UniqueDeviceId};
	% then remove one isntance of each
	[tmp i] = unique(useuid);
	useuid(i) = [];
	% and then don't use any of the remaining items (i.e. anything that was
	% present more than once)
	useuid = ~ismember({sws.UniqueDeviceId},useuid);
	for i = 1:length(sws)

		% Sometimes the "unique device ID" field isn't actually unique.
		% When it is, we use it, otherwise we use the section ID:
		if(useuid(i))
			sw(i).Name = sws(i).UniqueDeviceId;
		else
			sw(i).Name = sws(i).SectionId;
		end

		% object of the switch in this case is the line/section
		sw(i).SwitchedObj = ['Line.' sws(i).SectionId];

		if sws(i).SwitchIsOpen
			sw(i).Action = 'Open';
		else
			sw(i).Action = 'Close';
		end

		%TODO: check what kind of data type Lock field receives
		sw(i).Lock = logical(sws(i).SwitchIsTie);

		if(~sws(i).NearFromNode)
			sw(i).SwitchedTerm = 2;
		end
		
	end
	toc(t_);
	fns_handled{end+1} = 'InstSwitches';
end

%% Fuses
if(isfield(d,'InstFuses'));
	fprintf('Converting fuses... '); t_ = tic;
	fs = d.InstFuses;
	f(length(fs)) = dssfuse;

	for i = 1:length(fs)

		% Name without space
		f(i).Name = fs(i).UniqueDeviceId;

		if(fs(i).FuseIsOpen), f(i).Action = 'Open';
		else f(i).Action = 'Close';
		end

		% TODO: check if it can actually work when RatedCurrent is actual phase amps
		f(i).RatedCurrent = str2double(fs(i).AmpRating);

		% Handle spaces in sectionId
		f(i).MonitoredObj = ['Line.' fs(i).SectionId];
		if(~fs(i).NearFromNode)
			f(i).SwitchedTerm = 2;
		end

	end
	toc(t_);
	fns_handled{end+1} = 'InstFuses';
end

%% Transformers
if(isfield(d,'InstPrimaryTransformers'))
	fprintf('Converting transformers... '); t_ = tic;
	trs = d.InstPrimaryTransformers;
	tr(length(trs)) = dsstransformer;

	[busid, idx_bid] = ismember({trs.SectionId},{d.InstSection.SectionId});

	for i = 1:length(trs)

		tr(i).Name = trs(i).UniqueDeviceId;
		% find the matching line segment
		l_idx = find(strcmp({l.Name}',trs(i).SectionId));
		if(isempty(l_idx))
			tbus = regexprep(tr(i).Buses,'(\.\d)+$','');
			lbus = regexprep([{l.bus1}',{l.bus2}'],'(\.\d)+$','');
			l_idx = find(all([strcmp(tbus{1},lbus(:,1)) strcmp(tbus{2},lbus(:,2))],2));
		end
		tr(i).Phases = l(l_idx).Phases;

		% Buses
		tr(i).Buses = {d.InstSection(idx_bid(i)).FromNodeId d.InstSection(idx_bid(i)).ToNodeId};
		if(~trs(i).HighSideNearFromNode), tr(i).Buses = tr(i).Buses(end:-1:1); end
		%tr(i).Phases = trs(i).ConnectedPhases;

		% Set connections
		tr(i).Conns = {trs(i).HighSideConnectionCode ...
			trs(i).LowSideConnectionCode};

		% we only know the properties of the low winding, so we skip
		% directly to set those
		tr(i).wdg = 2;
		tr(i).Kv = trs(i).SpecNomKv;

		% Remove the line segment so that the transformer can actually do
		% something
		l(l_idx) = [];
		
	end
	toc(t_);
	fns_handled{end+1} = 'InstPrimaryTransformers';
end

%% Regulator Controllers (should be defined after transformers)
% ASSUMPTION: All raw transformer raw data is stored in 'trs' variable
if(isfield(d,'InstRegulators'))
	fprintf('Converting regulator controllers... '); t_ = tic;
	rs = d.InstRegulators;
	% If the database has a table describing regulator devices, issue a
	% warning to let the user know that our script doesn't work with those
	% at this time.
	if(isfield(d,'DevRegulators'))
		[haveRdevs haveRdevs] = ismember({rs.RegulatorType},{d.DevRegulators.RegulatorName});
	else
		haveRdevs = false(size(rs));
	end
	
	r(length(rs)) = dssregcontrol;
	
	[busid, idx_bid] = ismember({rs.SectionId},{d.InstSection.SectionId});
	
	if ~exist('trs','var')
		tr(length(rs)) = dsstransformer;
		tid = [];
		% initialize transformer index;
		tidex = 1;
	else
		[busid, tid] = ismember(trs.SectionId,{d.InstRegulators.SectionId});
		% number of existing transformer that match 
		ntr = sum(nonzeros(busid));
		
		tr( (length(tr)+1) : (length(tr)+ length(rs) -ntr) ) = dsstransformer;
		
		% initialize transformer index;
		tidex = length(tr) + 1;
	end
	
	for i = 1:length(rs)
		ri = r(i);
		rsi = rs(i);
		
		% create transformer if not exists
		if ~ismember(i,tid)
			tri = tr(tidex);
			% Assign basic properties for a two-winding transformer
			tri.Name = rsi.UniqueDeviceId;
			tri.Buses = {d.InstSection(idx_bid(i)).FromNodeId d.InstSection(idx_bid(i)).ToNodeId};
			tri.Phases = rsi.ConnectedPhases;
			% don't forget to specify connection style for both windings
			tri.Conns = {rsi.RegulatorConfiguration, rsi.RegulatorConfiguration};
			% select the winding that we're going to put some taps on
			if(rsi.TapsNearFromNode == 1)
				tri.wdg = 1;
			else
				tri.wdg = 2;
			end
			% Fill in some additional transformer properties if they're
			% present in the model
			if(haveRdevs(i))
				RegDev = d.DevRegulators(haveRdevs(i));
				tri.MaxTap = 1+RegDev.RaiseAndLowerMaxPercentage/100;
				tri.MinTap = 2-tri.MaxTap;
				tri.NumTaps = RegDev.NumberOfTaps;
				tri.kVAs = RegDev.RegulatorRatedKva*[1 1];
				tri.kVs = RegDev.RegulatorRatedVoltage*[1 1];
				% I don't think the following two are quite right
				warning('dssconversion:attention','Using some regulator settings that haven''t been tested');
				tri.Rs = RegDev.PercentZOnRegulatorBase*[1 1];
				tri.imag = tri.R *RegDev.RegulatorXRRatio;
			else
				% Mostly the defaults are at least as good as what I can do
				% because I'd be making it up, but there are a couple:
				% if the tap limiter is being used, its range is a
				% reasonable number of taps to use
				% This also doesn't let us guess maxtaps or mintaps because
				% they're pu, not integer.
				if(rsi.TapLimiterActive)
					tri.NumTaps = rsi.TapLimiterHighSetting - rsi.TapLimiterLowSetting;
				end
				% we can try to extract the current/kv ratings from the
				% regulator type name.  In the examples I've seen, these
				% are often supplied as, e.g. '200A 12' or '12.47KV 350'.
				% But this will just end up being a guess
				kva = regexp(rsi.RegulatorType,'([\d.]+)A ([\d.]+)','tokens','once');
				if(isempty(kva))
					kva = regexpi(rsi.RegulatorType,'([\d.]+)KV ([\d.]+)','tokens','once');
					kva = kva(end:-1:1);
				end
				if(~isempty(kva))
					tri.kVAs = prod(str2double(kva))*[1 1];
					tri.kVs = str2double(kva{2})*[1 1];
				else
					% and finally, if all that's failed, try to use the
					% feeder voltage
					try
						tri.kVs = d.InstFeeders.NominalKvll*[1 1];
						if(strcmp(tri.Conn,'wye'))
							tri.kVs = tri.kVs/sqrt(3);
						end
					catch
					end
				end
			end
			% grab the phase info; 'phase' var is used further down, so
			% don't mess with it.
			if(rsi.TapsAreGangOperated)
				phase = rsi.GangOperatingPhase;
				% let the transformer object convert the integer tap
				% position into a pu tap position; it does this when we use
				% integer data types
% 				tri.tap = int16(rsi.(['TapPositionPhase' phasemap.(phase)(2)]));
			else
				phase = rsi.ConnectedPhases;
				if(sum(phase~=' ' & phase~='N')>1), warning('dssconversion:ThreePhase','Converting a multiphase regulator that is not gang operated may lose info as implemented'); end
				% adopt the most common tap position as the one to use
				tri.tap = int16(mode(nonzeros([rsi.TapPositionPhase1 rsi.TapPositionPhase2 rsi.TapPositionPhase3])));
			end

			tr(tidex) = tri;
			tidex = tidex + 1;
		else
			% load old tranformer
			tri = tr(ismember(rsi.SectionId,{trs.SectionId}));
			% This really shouldn't ever happen in my understanding of the
			% SynerGee model format.  I'm therefore opting not to spend
			% time writing the code correctly to fill out the transformer
			% in this case.
			warning('dssconversion:duplicateRegulator','Found an existing transformer with the same name ''%s'' used by a regulator.',tri.Name);
		end
		
		% Start setting properties on the regcontrol object
		ri.Name = rsi.UniqueDeviceId;
		ri.transformer = tri.Name;
		
		% Get bus ID from SectionID ('RG 05201317_05201317A') and
		% NearFromNode properties
		if rsi.NearFromNode == 1
			bus = d.InstSection(idx_bid(i)).FromNodeId;
		else
			bus = d.InstSection(idx_bid(i)).ToNodeId;
		end
		% select the Phases to monitor on that bus
		% doing it this way actually implies something like Synergee's
		% TapsAreGangOperated, so we'd want to look back at that if they're
		% not.
		[x y phase] = dataclean(phase,'monitoredphase');
		ri.PTphase = phase;
		%bus = [bus '.' phase];
		%ri.bus = bus;
		
		% select the appropriate winding for the regulator controller
		if rsi.TapsNearFromNode == 1
			ri.winding = 1;
		else
			ri.winding = 2;
		end
		
		% Set (or guess) the pt ratio
		if(haveRdevs(i))
			ri.ptratio = RegDev.PTRatio;
			% I think these two CT numbers are the same thing, but OpenDSS
			% doesn't document very well what it's doing.  Based on reading
			% the code, I think it uses this as the ratio between line
			% current and control current, which is exactly what synergee
			% does
			ri.CTprim = RegDev.CTRating;
			
		else
			% PT ratio is supposed to be assigned so that the nominal
			% line-line voltage corresponds to a control Line-Neutral
			% voltage of 120 afaict
			ri.ptratio = tri.kVs(1)/sqrt(3)/120;
			% skip CT rating I guess?  I'm not sure how I'd pick one, and
			% it shouldn't matter anyway as long as we don't use the LDC
			% features in OpenDSS.
		end
		
		% Use either the specified control values (for gang operated taps)
		% or the most common control values for the forward direction
		if(rsi.TapsAreGangOperated)
			phase = phasemap.(rsi.GangOperatingPhase)(2);
			ri.band = rsi.(['ForwardBWDialPhase' phase]);
			ri.vreg = rsi.(['ForwardVoltageSettingPhase' phase]);
		else
			ri.band = mode(nonzeros([rsi.ForwardBWDialPhase1 rsi.ForwardBWDialPhase2 rsi.ForwardBWDialPhase3]));
			ri.vreg = mode(nonzeros([rsi.ForwardVoltageSettingPhase1 rsi.ForwardVoltageSettingPhase2 rsi.ForwardVoltageSettingPhase3]));
		end
		if(any([rsi.ForwardRDialPhase1 rsi.ForwardRDialPhase2 rsi.ForwardRDialPhase3]))
			warning('dssconversion:unimplemented','Regulator %s uses LDC, which is not implemented',ri.Name);
		end
		
		% convert the reverse control mode
		% see chart on p167 of synergee tech manual
		switch(rsi.ReverseSensingMode)
			case 'LF'%Lock Forward
			case 'IN'%Neutral Idle
				ri.revNeutral = 'Yes';
			case {'NR','CG'} %no reverse or cogeneration?
			case {'IR','BD','LR'}% idle reverse
				warning('dssconversion:unimplemented','I don''t do anything with the reverse sensing mode %s yet',rsi.ReverseSensingMode)
			otherwise % means I guessed the acronyms above incorrectly
				error('dssconversion:implfailure','Bad Reverse Sensing mode: %s',rsi.ReverseSensingMode);
		end
		% setting for reverse direction operation
		if any([rsi.ReverseVoltageSettingPhase1 rsi.ReverseVoltageSettingPhase2 rsi.ReverseVoltageSettingPhase3])
			ri.reversible = 'Yes';
			% grab the voltage band
			if(rsi.TapsAreGangOperated)
				ri.revband = rsi.(['ReverseBWDialPhase' phase]);
				ri.revvreg = rsi.(['ReverseVoltageSettingPhase' phase]);
			else
				ri.revband = mode(nonzeros([rsi.ReverseBWDialPhase1 rsi.ReverseBWDialPhase2 rsi.ReverseBWDialPhase3]));
				ri.revvreg = mode(nonzeros([rsi.ReverseVoltageSettingPhase1 rsi.ReverseVoltageSettingPhase2 rsi.ReverseVoltageSettingPhase3]));
			end
		end
		
		if(~(rsi.RegulatorIsOn) || rsi.ManualTapOperation)
			ri.maxtapchange = 0; %disable tap changes
		end
		if(rsi.FirstHouseActive)
			ri.vlimit = rsi.FirstHouseHighVoltage;
			if(rsi.FirstHouseLowVoltage~=0)
				warning('dssconversion:fhlv','First House Low Voltage protection does not exist in opendss');
			end
		end
		
		if(rsi.NomKvMult~=1)
			warning('dssconversion:unimplemented','Regulator uses a nonunity kv multiplier');
		end
		
		r(i) = ri;		
	end
	
	%TODO: finally, remove lines that represent these regulator(s)
	[val l_idx] = ismember(dataclean({d.InstRegulators.SectionId},'name'),{l.Name});
	l(nonzeros(l_idx)) = [];

	toc(t_);
	fns_handled{end+1} = 'InstRegulators';
end

%% Reclosers: 
% TODO: Need for fault study (need more careful research accurate
% implemetation of reclosers)
%
% Some parameters that we are not sure how to implement correctly are:
% - 'ampRating' (630)
% - FastPhaseCurve: '101'
% - SlowPhaseCurve: '101'
% - FastGroundCurve: '101'
% - SlowGroundCurve: '101'
% - InterruptRatingAmps: 12500
if(isfield(d,'InstReclosers'))
	fprintf('Converting reclosers... '); t_ = tic;
	rcs = d.InstReclosers;
	rc(length(rcs)) = dssrecloser;

	for i = 1:length(rcs)

		rc(i).Name = rcs(i).UniqueDeviceId;
		rc(i).MonitoredObj = ['Line.' rcs(i).SectionId];
		if(~rcs(i).NearFromNode)
			rc(i).SwitchedTerm = 2;
		end
		rc(i).NumFast = rcs(i).FastPhaseCount;
		rc(i).Shots = rcs(i).FastPhaseCount + rcs(i).SlowPhaseCount;
		rc(i).TDGrDelayed = rcs(i).SlowGroundTimeAddSec;
		rc(i).TDGrFast = rcs(i).FastGroundTimeAddSec;
		rc(i).TDPhDelayed = rcs(i).SlowPhaseTimeAddSec;
		rc(i).TDPhFast = rcs(i).FastPhaseTimeAddSec;
	end
	toc(t_);
	fns_handled{end+1} = 'InstReclosers';
end

%% Circuit/Feeder
% According to the SynerGEE documentation, there are two kinds of sources
% in synergee.  A "feeder" is used when substation data is unavailable, and
% is "usually setup to represent some point near the secondary of the
% substation transformer."  A "substation" had historically included the
% transformer as well, but that data has been moved to the primary
% transfomers table when applicable, leaving the "substation" very similar
% to the feeder.  Not sure how "substation" is stored in a table since we
% don't have an example of one.
cir = dsscircuit;
fi = d.InstFeeders;
cir.basekv = fi.NominalKvll;
cir.Name = fi.SubstationId;
cir.bus1 = fi.FeederId;
cir.R1 = fi.PosSequenceResistance;
cir.X1 = fi.PosSequenceReactance;
cir.R0 = fi.ZeroSequenceResistance;
cir.X0 = fi.ZeroSequenceReactance;
if(mod(fi.ByPhVoltDegPh1 - fi.ByPhVoltDegPh2,360) > 0) % Phase 1 leads
	cir.Sequence = 'pos';
else
	cir.Sequence = 'neg';
end
fns_handled{end+1} = 'InstFeeders';

%% Bus list with locations
if(isfield(d,'Node'))
	fprintf('Generating buslist... '); t_ = tic;

	output.buslist.id = dataclean({d.Node.NodeId},'name')';
	output.buslist.coord = [vertcat(d.Node.X), vertcat(d.Node.Y)];
	fns_handled{end+1} = 'Node';

	toc(t_);
end

%% Check what fields are left
fns_handled = [fns_handled {'InstMymLargeCust','InstMymLoads','InstMymMeter','InstViews','SAI_Control'}]; % to specify others all at once
fn_unknown = setdiff(fieldnames(d),fns_handled);
disp('Unknown object types:');
disp(fn_unknown);

%% Output
voltages = [];
if(exist('cir','var'))
	output.circuit = cir;
	voltages = [voltages cir.basekv];
end
if(exist('capcon','var'))
	output.capcontrol = capcon;
end
if(exist('ld','var'))
	output.load = ld;
	voltages = [voltages ld.Kv];
end
if(exist('ge','var'))
	output.generator = ge;
end
if(exist('l','var'))
	output.line = l;
end
if(exist('lc','var'))
	output.linecode = lc;
end
if(exist('sw','var'))
	output.switch = sw;
end
if(exist('tr','var'))
	output.transformer = tr;
	voltages = [voltages tr.kVs];
end
if(exist('r','var'))
	output.regcontrol = r;
end
if(exist('f','var'))
	output.fuse = f;
end
if(exist('rc','var'))
	output.recloser = rc;
end
if(exist('cp','var'))
	output.capacitor = cp;
	voltages = [voltages cp.Kv];
end
if(~isempty(voltages))
	output.basevoltages = sort(unique(voltages),'descend');
end
end
