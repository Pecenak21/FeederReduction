function obj = dssrecloser(s,varargin)
% OpenDSS recloser object
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
	obj.defaults = struct('Name',{''}, ...
		'Like', {''}, ...
		'Action','', ... %{Trip/Open | Close}  Action that overrides the Recloser control. Simulates manual control on recloser "Trip" or "Open" causes the controlled element to open and lock out. "Close" causes the controlled element to close and the Recloser to reset to its first operation.
		'Delay',0, ...
		'enabled','Yes', ... %{Yes|No or True|False} Indicates whether this element is enabled.
		'GroundDelayed','none', ... %Name of the TCC Curve object that determines the Ground Delayed trip.  Must have been previously defined as a TCC_Curve object. Default is none (ignored).Multiplying the current values in the curve by the "groundtrip" value gives the actual current.
		'GroundFast','none', ... %Name of the TCC Curve object that determines the Ground Fast trip.  Must have been previously defined as a TCC_Curve object. Default is none (ignored). Multiplying the current values in the curve by the "groundtrip" value gives the actual current.
		'GroundInst','none', ... %Actual amps for instantaneous ground trip which is assumed to happen in 0.01 sec + Delay Time.Default is 0.0, which signifies no inst trip.
		'GroundTrip',0, ... %Multiplier or actual ground amps (3I0) for the ground TCC curve.  Defaults to 1.0.
		'MonitoredObj','', ...
		'MonitoredTerm',1, ... %Number of the terminal of the circuit element to which the Recloser is connected. 1 or 2, typically.  Default is 1.
		'NumFast',1, ... %Number of Fast (fuse saving) operations.  Default is 1. (See "Shots")
		'PhaseDelayed','D', ... %Name of the TCC Curve object that determines the Phase Delayed trip.  Must have been previously defined as a TCC_Curve object. Default is "D".Multiplying the current values in the curve by the "phasetrip" value gives the actual current.
		'PhaseFast','A', ... %Name of the TCC Curve object that determines the Phase Fast trip.  Must have been previously defined as a TCC_Curve object. Default is "A". Multiplying the current values in the curve by the "phasetrip" value gives the actual current.
		'PhaseInst',0, ... %Actual amps for instantaneous phase trip which is assumed to happen in 0.01 sec + Delay Time. Default is 0.0, which signifies no inst trip. 
		'PhaseTrip',1, ... %Multiplier or actual phase amps for the phase TCC curve.  Defaults to 1.0.
		'RecloseIntervals',[0.5 2 2], ... %Array of reclose intervals.  Default for Recloser is (0.5, 2.0, 2.0) seconds. A locked out Recloser must be closed manually (action=close). 
		'Reset',15, ... %Reset time in sec for Recloser.  Default is 15. 
		'Shots',4, ... %Total Number of fast and delayed shots to lockout.  Default is 4. This is one more than the number of reclose intervals.
		'SwitchedObj','',... %Name of circuit element switch that the Recloser controls. Specify the full object name.Defaults to the same as the Monitored element. This is the "controlled" element.
		'SwitchedTerm',1,... %Number of the terminal of the controlled element in which the switch is controlled by the Recloser. 1 or 2, typically.  Default is 1.
		'TDGrDelayed',1,... %Time dial for Ground Delayed trip curve. Multiplier on time axis of specified curve. Default=1.0.
		'TDGrFast',1,... %Time dial for Ground Fast trip curve. Multiplier on time axis of specified curve. Default=1.0.
		'TDPhDelayed',1,... %Time dial for Phase Delayed trip curve. Multiplier on time axis of specified curve. Default=1.0.
		'TDPhFast',1,... %Time dial for Phase Fast trip curve. Multiplier on time axis of specified curve. Default=1.0.
		'Basefreq',60);
	obj.fieldnames = fieldnames(obj.defaults);
	% set the data structure to be a blank version of the defaults
	% structure
	obj.data = obj.defaults;
	obj.data(1) = []; obj.data(1).Name = ''; obj.data(1).Like = '';
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