function obj = dssenergymeter(s,varargin)
% OpenDSS energy meter object
% This class keeps track of default values for each of the parameters of
% the opendss object, as well as user-assigned parameters, and then can
% pass them on to other objects for use in converting opendss files to
% other formats (see get() ) or can print itself as a char string (see
% char() ) for converting to OpenDSS from other formats.
%
% Generally, for converting to opendss you will want to create one or more
% empty objects and then fill in some values, then export to a text file
% using the char() function.
%
% There are three valid ways to create these objects (specified generically
% so replace 'dssobject' with the class name):
%	dssobject()
%	dssobject(struct)
%	dssobject('prop1',value1,'prop2',value2,...)
%
% The first creates an empty object ready to be filled in
% The second and third populate as many values as possible based on the
% data passed in as arguments or as a struct.  Unknown properties are
% ignored.  Order of properties is not important in the third format.
% Unfortunately, we currently do not support struct arrays or the cell
% array style input that struct() accepts for creating struct arrays.

% The class structure is as follows:
% obj.defaults contains the default values of all the parameters
% obj.data contains user-set values of all parameters
% obj.fieldnames contains the fieldnames of the two structs to save you
% making lots of extra calls.
%
% The concept of this separation of data is that it allows us to present a
% "complete, as opendss would use it" version of the model for converting
% to other formats while only setting the fields we need to explicitly when
% converting TO OpenDSS.
% 
% The question has come up about what to do if there are values we'd like
% to keep in a general conversion utility that don't neatly fit into
% openDSS.  First, make sure they actually don't (openDSS has a lot of
% stuff, and we keep finding ways that it actually does support the
% property we thought it didn't).  The idea I had was to add a third
% sub-struct to hold custom data in some format.  Then modify the 'get' and
% 'set' accessors to use data from it when appropriate.  That way, we can
% continue to maintain separation of data.
%
% When duplicating this class, it should only be necessary to change a few
% methods.  You'll want to change the object defaults below in this file,
% and you'll want to change how the object deals with handling new data in
% the set() function, but that should be about it.  (of course, as you add
% new things, you may also find bugs in various other files that need
% fixing).

% get the class name, for use in various places
cn = mfilename('class');
% if we were asked to create a copy of an existing object of our class, do
% that now
if(nargin==1 && isa(s,cn))
	obj = struct(s);
	done = 1;
else %otherwise assign default values
	obj.defaults = struct('Name','', ...
		'Like', '', ...
		'threephaseLosses','YES',... %{Yes | No}  Default is YES. Compute Line losses and segregate by 3-phase and other (1- and 2-phase) line losses. 
		'action','',... 
			%{
			Clear (reset) | Save | Take | Zonedump | Allocate | Reduce 
			(A)llocate = Allocate loads on the meter zone to match PeakCurrent.
			(C)lear = reset all registers to zero
			(R)educe = reduces zone by merging lines (see Set Keeplist & ReduceOption)
			(S)ave = saves the current register values to a file.
			   File name is "MTR_metername.CSV".
			(T)ake = Takes a sample at present solution
			(Z)onedump = Dump names of elements in meter zone to a file
			   File name is "Zone_metername.CSV".
			%}
		'element','',...%Name (Full Object name) of element to which the monitor is connected.
		'enabled','',... %{Yes|No or True|False} Indicates whether this element is enabled.
		'kVAemerg',0,...%Upper limit on kVA load in the zone, Emergency configuration. Default is 0.0 (ignored). Overrides limits on individual lines for overload UE. With "LocalOnly=Yes" option, uses only load in metered branch.
		'kVAnormal',0,...%Upper limit on kVA load in the zone, Normal configuration. Default is 0.0 (ignored). Overrides limits on individual lines for overload EEN. With "LocalOnly=Yes" option, uses only load in metered branch.
		'LineLosses','YES',...%{Yes | No}  Default is YES. Compute Line losses. If NO, then none of the losses are computed.
		'LocalOnly','NO',... %{Yes | No}  Default is NO.  If Yes, meter considers only the monitored element for EEN and UE calcs.  Uses whole zone for losses.
		'Losses','YES',... %{Yes | No}  Default is YES. Compute Zone losses. If NO, then no losses at all are computed.
		'Mask',1,... % Mask for adding registers whenever all meters are totalized.  Array of floating point numbers representing the multiplier to be used for summing each register from this meter. Default = (1, 1, 1, 1, ... ).  You only have to enter as many as are changed (positional). Useful when two meters monitor same energy, etc.
		'Option','E',...
			%{
				Enter a string ARRAY of any combination of the following. Options processed left-to-right:

				(E)xcess : (default) UE/EEN is estimate of energy over capacity 
				(T)otal : UE/EEN is total energy after capacity exceeded
				(R)adial : (default) Treats zone as a radial circuit
				(M)esh : Treats zone as meshed network (not radial).
				(C)ombined : (default) Load UE/EEN computed from combination of overload and undervoltage.
				(V)oltage : Load UE/EEN computed based on voltage only.

				Example: option=(E, R)
			%}
		'PeakCurrent',[400 400 400],... %ARRAY of current magnitudes representing the peak currents measured at this location for the load allocation function.  Default is (400, 400, 400). Enter one current for each phase
		'PhaseVoltageReport','NO',... %{Yes | No}  Default is NO.  Report min, max, and average phase voltages for the zone and tabulate by voltage base. Demand Intervals must be turned on (Set Demand=true) and voltage bases must be defined for this property to take effect. Result is in a separate report file.
		'SeqLosses','YES',... %{Yes | No}  Default is YES. Compute Sequence losses in lines and segregate by line mode losses and zero mode losses.
		'terminal',[],... % Number of the terminal of the circuit element to which the monitor is connected. 1 or 2, typically.
		'VbaseLosses','YES',... %{Yes | No}  Default is YES. Compute losses and segregate by voltage base. If NO, then voltage-based tabulation is not reported.
		'XfmrLosses','YES',... %{Yes | No}  Default is YES. Compute Transformer losses. If NO, transformers are ignored in loss calculations.
		'Zonelist',{''},...
			%{
				ARRAY of full element names for this meter's zone.  Default is for meter to find it's own zone. If specified, DSS uses this list instead.  Can access the names in a single-column text file.  Examples: 
				zonelist=[line.L1, transformer.T1, Line.L3] 
				zonelist=(file=branchlist.txt)
			%}
		'Basefreq',60 ...
		);
	obj.fieldnames = fieldnames(obj.defaults);
	% set the data structure to be a blank version of the defaults
	% structure
	obj.data = obj.defaults;
	obj.data(1) = []; obj.data(1).Name = ''; obj.data(1).Like = '';
	% some fields require re-naming
	obj.namemap = struct('threephaseLosses','3phaseLosses');
end

% make it an object
obj = class(obj,cn);

% if we had struct or variable input arguments, now's the time to use them.
% just loop through all the values and call set() with them.
if(nargin == 0 || exist('done','var'))
	return;
elseif(nargin == 1 && isstruct(s))
	% here we extract a list of fieldnames and corresponding values to use
	% for setting properties
	fn = fieldnamefix(fieldnames(s),obj.fieldnames);
	v = struct2cell(s);
elseif(mod(nargin,2)==0)
	% Similarly if inputs are passed as key-value pairs, we separate them
	% by keys and values
	v = reshape([{s},varargin(:)'],[2 nargin/2])';
	fn = fieldnamefix(v(:,1),obj.fieldnames);
	v = v(:,2);
else
	error('Invalid inputs');
end
% loop over the keys/values to set them
for i=1:length(fn)
	if(isempty(fn{i})), continue; end
	obj = set(obj,fn{i},v{i});
end

end