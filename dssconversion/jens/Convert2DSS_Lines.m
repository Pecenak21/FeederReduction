function l = Convert2DSS_Lines(s,lc)
% l = Convert2DSS_Lines(s)
%
% PURPOSE : Converts lines from object to OpenDSS
%
%
% INPUT :   s: struct that has all entries (e.g., o.R1)
%
% OUTPUT :  OpenDSS Lines object
%           see page 106 in OpenDSS manual V7.4.3(March 2012)

% l(length(s)) = dssline;
% make phase mapping easier:
if isfield(s,'Phase1')
    phasemap.(s.Phase1) = '.1';
else
    phasemap.A = '.1';
end
if isfield(s,'Phase2')
    phasemap.(s.Phase2) = '.2';
else
    phasemap.B = '.2';
end
if isfield(s,'Phase2')
    phasemap.(s.Phase3) = '.3';
else
    phasemap.C = '.3';
end
phasemap.N = '.0';
    


lineunits = lc.Units;
% if(any(strcmpi(d.SAI_Control.LengthUnits,{'English','English2'})))
% 	lineunits = 'Ft';
% elseif(strcmpi(d.SAI_Control.LengthUnits,'Metric'))
% 	lineunits = 'meter';
% end

missing_linecodes = {};
for i = 1:length( s )
    l(i)=dssline;
	l(i).Name = s(i).SectionID;
	% setting the linecode
    % 1. search for the linecode id
    
	%if strcmp(s(i).SectionType,'other')~=1 % only process UG or OH line types (no switches)
	%[x idx] = ismember(dataclean(s(i).LineCableID,'name'), get(lc(:),'Name'));
	if strcmp(s(i).SectionType,'UG') || strcmp(s(i).SectionType,'OH')
        
		idx=GetIndex_StrInCell(lc(:).Name,dataclean(s(i).LineCableID,'name'));
		if isempty(idx)
			% check for missing
			if(~isempty(s(i).LineCableID))
				missing_linecodes = [missing_linecodes, s(i).LineCableID];
            end
			% fill in using the default
			idx=GetIndex_StrInCell(lc(:).Name,'DEFAULT');
        end
        if length(idx)==2
			if strcmp(s(i).SectionType,'UG') % two default values, one for OH (first) and one for UG (second)
				idx = idx(2);
			else
				idx = idx(1);
			end
		end
		% 2. set linecode of the line to be that linecode object
		l(i).LineCode = lc(idx);
        % do length
		l(i).Length = s(i).Length;
    elseif strcmp(s(i).SectionType,'geometry')
        l(i).geometry=[s.SpacingID '__' s.CondID_A '__' s.CondID_B '_' s.CondID_C '__' s.CondID_N];
        l(i).Length = s(i).Length;
    elseif strcmp(s(i).SectionType,'switch') || strcmp(s(i).SectionType,'other')|| strcmp(s(i).SectionType,'capacitor')
        l(i).switch='y';
        if ~strcmp(s(i).ClosedPhase,'NONE')
            l(i).enabled='true';
        else
            l(i).enabled='false';
        end
        l(i).Length = 0.001;
	end

	%only report a neutral phase if we're going to ground it
	% 	secphases = s(i).SectionPhases;
	secphases = s(i).Phase;
	secphases(secphases==' ') = [];
% js: commented out this if statement, i.e., never report a neutral, may
% have to change this later if neutral connections are important
% 	if(~s(i).NeutIsGrounded)
		secphases(secphases=='N') = [];
% 	end
	l(i).Phases = length(secphases);
	
	secphstr = '';
	for phase_idx = 1:length(secphases)
		secphstr = [secphstr phasemap.(secphases(phase_idx))];
	end
	l(i).bus1 = [s(i).FromNodeID secphstr];
	l(i).bus2 = [s(i).ToNodeID secphstr];
	%if(~s(i).IdenticalPhaseConductors)
	%	warning('dssconversion:unimplemented','Use of different phases for different conductors is not implemented');
	%end
	l(i).Units = lineunits;
	%end
end
if(~isempty(missing_linecodes))
	missing_linecodes = unique(missing_linecodes);
	warning('dssconversion_CYME:missing_linecode','Could not find the following linecodes: \n%s', sprintf('%s\n',missing_linecodes{:}));
end


end
% obj.defaults = struct( ...
% 		'Name','',...
% 		'Geometry','',... % Geometry code for LineGeometry Object. Supercedes any previous definition of line
%                           impedance. Line constants are computed for each frequency change or rho change.
%                           CAUTION: may cause the number of phases to be redefined.
% 		'LineCode','',... % Name of an existing LineCode object containing impedance definitions.
% 		'Phases',3,...
% 		'bus1','',...
% 		'bus2','',...
% 		'Units','None',...	% Units can be used to allow specification of impedance and line length in different units.  See below.
% 		... % we want to make sure that units are printed before all of the quantities they apply to, so that the values won't accidentally get converted to a new system of units when we read in the file
% 		'Length',1,... %kft, Length multiplier to be applied to the impedance data.
% 		'Switch','no',... % switch effects all the R properties, so we place it before them so that if they are set too, they will be written after and not get lost
% 		'R1',0.058,... %ohms per 1000 ft
% 		'X1',0.1206,... %ohms per 1000 ft
% 		'R0',0.1784,... %ohms per 1000 ft
% 		'X0',0.4047,... %ohms per 1000 ft
% 		'C1',3.4,... %nF per 1000 ft
% 		'C0',1.6,... %nF per 1000 ft
% 		'BaseFreq'  ,60     ,...
% 		'Normamps'  ,[]     ,...
% 		'Emergamps' ,[]     ,...
% 		'Faultrate' ,0.0005     ,...
% 		'Pctperm'   ,[]     ,...
% 		'Rg'        ,0.01805     ,...
% 		'Xg'        ,0.155081     ,...
% 		'Rho'       ,100     ,...
% 		'Like'      ,''     ,...
% 		'Repair',[],... %hours to repair
% 		'Rmatrix',[],...
% 		'Xmatrix',[],...
% 		'Cmatrix',[],...
% 		'EarthModel','',...
% 		'enabled','' ...
% 		);
	% Units: The default assumption ('None') for units implies that the
	% impedances are specified per unit length in the same units as the
	% line length.  OpenDSS allows you to specify any of a variety of units
	% (miles, feet, kfeet, cm, m, km), which are matched based on the first
	% two letters (or one letter for 'm').  Unit conversion only works well
	% in the case where you are applying a linecode with one set of units
	% (e.g. kft; so impedances in ohms/kft) to a line with different units
	% (e.g. ft).  In this case, the length is converted to the same units
	% as impedance when calculating the Y matrix, or the impedances are
	% converted to the same units as the length (and the 'units' param of
	% the line) for display (by OpenDSS; we're not that sophisticated here
	% yet).  Attempting to specify impedance and line length with different
	% units for just one line object is tricky.  Conversions are applied as
	% before (i.e. to the length for calculation purposes or the impedance
	% for display purposes) however _setting_ the unit is complicated:
	%	* The first time a unit is set serves only to adjust the 'current'
	%	unit, and does not effect the scale factor, because opendss
	%	refuses to convert from 'none' units to anything else
	%	* setting any of the impedance values or a linecode resets the
	%	units to none
	%	* by applying the units property a SECOND time, you can achieve the
	%	same effect as specifying a linecode with the first set of units
	%	and a line with the second set.  e.g. to interpret impedances in
	%	ohms/kft and lengths in ft, one might do:
	%	New line.myline R=0.05 length=534 units=kft units=ft
	%	Unfortunately, the current matlab class architecture we've setup is
	%	incapable of specifying the same parameter twice, so for creating
	%	OpenDSS files that behave like this, you either need to use
	%	linecodes or else just convert the units yourself.
