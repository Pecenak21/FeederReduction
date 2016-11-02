function obj = dssregcontrol(s,varargin)
% OpenDSS regulator control object
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
	obj.defaults = struct(...
		'Name',{''}, ...
		'Like','',...
		'transformer','', ... %Name of Transformer element to which the RegControl is connected. Do not specify the full object name; "Transformer" is assumed for the object class. Example: Transformer=Xfmr1
		'winding',1, ... %Number of the winding of the transformer element that the RegControl is monitoring. 1 or 2, typically. Side Effect: Sets TAPWINDING property to the same winding.
		'vreg',120, ... %Voltage regulator setting, in VOLTS, for the winding being controlled. Multiplying this value times the ptratio should yield the voltage across the WINDING of the controlled transformer. Default is 120.0
		'band',3.0, ... %Bandwidth in VOLTS for the controlled bus (see help for ptratio property). Default is 3.0
		'delay',15, ... %Time delay, in seconds, from when the voltage goes out of band to when the tap changing begins. This is used to determine which regulator control will act first. Default is 15. You may specify any floating point number to achieve a model of whatever condition is necessary.
		'ptratio',60, ... %Ratio of the PT that converts the controlled winding voltage to the regulator voltage. Default is 60. If the winding is Wye, the line?to?neutral voltage is used. Else, the line?to?line voltage is used.
		'CTprim',0.2, ... %Rating, in Amperes, of the primary CT rating for converting the line amps to control amps.The typical default secondary ampere rating is 0.2 Amps (check with manufacturer specs).
		'R',0, ... % R setting on the line drop compensator in the regulator, expressed in VOLTS.
		'X',0, ... % X setting on the line drop compensator in the regulator, expressed in VOLTS.
		'PTPhase',1, ... %phase to monitor for voltage control
		'tapwinding',1, ... % Winding containing the actual taps, if different than the WINDING property. Defaults to the same winding as specified by the WINDING property.
		'bus','', ... % Name of a bus (busname.nodename) in the system to use as the controlled bus instead of the bus to which the transformer winding is connected or the R and X line drop compensator settings. Do not specify this value if you wish to use the line drop compensator settings. Default is null string. Assumes the base voltage for this bus is the same as the transformer winding base specified above. Note: This bus (1? phase) WILL BE CREATED by the regulator control upon SOLVE if not defined by some other device. You can specify the node of the bus you wish to sample (defaults to 1). If specified, the RegControl is redefined as a 1?phase device since only one voltage is used.
		'debugtrace','No', ... % {Yes | No* } Default is no. Turn this on to capture the progress of the regulator model for each control iteration. Creates a separate file for each RegControl named "REG_name.CSV".
		'EventLog','Yes', ...% {Yes/True* | No/False} Default is YES for regulator control. Log control actions to Eventlog.
		'inversetime','No', ...% {Yes | No* } Default is no. The time delay is adjusted inversely proportional to the amount the voltage is outside the band down to 10%.
		'maxtapchange',16,... % Maximum allowable tap change per control iteration in STATIC control mode. Default is 16. Set this to 1 to better approximate actual control action. Set this to 0 to fix the tap in the current position.
		'revband',[],... %Bandwidth for operating in the reverse direction.
		'revDelay',60,... %Time Delay in seconds (s) for executing the reversing action once the threshold for reversing has been exceeded. Default is 60 s.
		'reversible','No',... %{Yes |No*} Indicates whether or not the regulator can be switched to regulate in the reverse direction. Default is No.Typically applies only to line regulators and not to LTC on a substation transformer.
		'revNeutral','No',... %{Yes | No*} Default is no. Set this to Yes if you want the regulator to go to neutral in the reverse direction.
		'revR',[],... %R line drop compensator setting for reverse direction.
		'revThreshold',100,... % kW reverse power threshold for reversing the direction of the regulator. Default is 100.0 kw.
		'revvreg',[],... % Voltage setting in volts for operation in the reverse direction.
		'revX',[],... %X line drop compensator setting for reverse direction.
		'tapdelay',2,... %Delay in sec between tap changes. Default is 2. This is how long it takes between changes after the first change.
		'TapNum',0,... %An integer number indicating the tap position that the controlled transformer winding tap position is currently at, or is being set to.  If being set, and the value is outside the range of the transformer min or max tap, then set to the min or max tap position as appropriate. Default is 0
        'vlimit',0.0,... %Voltage Limit for bus to which regulated winding is connected (e.g. first customer). Default is 0.0. Set to a value greater then zero to activate this function.
        'enabled',''...
		);
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