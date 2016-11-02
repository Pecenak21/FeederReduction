function obj = dssvsource(s,varargin)
% OpenDSS circuit object
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

% from page 105 in OpenDSS manual V7.6(November 2012):
% Voltage source. This is a special power conversion element. It is special because voltage
% sources are used to initialize the power flow solution with all other injection sources set to zero.
% A Vsource object is a two?terminal, multi?phase Thevenin equivalent. That is, it is a voltage
% source behind an impedance. The data are specified as it would commonly be for a power
% system source equivalent: Line?line voltage (kV) and short circuit MVA.
% The most common way to use a voltage source object is with the first terminal connected to the
% bus of interest with the second terminal connected to ground (voltage reference). In this usage,
% the connection of the second terminal may be omitted. In 2009, the voltage source was changed
% from a single?terminal device to a two?terminal device. This allows for the connection

% get the class name, for use in various places
cn = mfilename('class');
% if we were asked to create a copy of an existing object of our class, do
% that now
if(nargin==1 && isa(s,cn))
	obj = struct(s);
	done = 1;
else %otherwise assign default values
	obj.defaults = struct(...
		'Name','SOURCE', ...
		'Like','', ...
		'Phases', 3, ...
		'bus1','', ... % Name of bus to which the source's first terminal is connected. Remember to specify the node order if the terminals are connected in some unusual manner. Side effect: The processing of this property results in the setting of the Bus2 property so that all conductors in terminal 2 are connected to ground.For example, Bus1= busname Has the side effect of setting Bus2=busname.0.0.0
		'bus2','', ... & Name of bus to which the source’s second terminal is connected. If omitted, the second terminal is connected to ground (node 0) at the bus designated by the Bus1 property.
        'basekv',12.47,... % base or rated Line?to?line kV.
		'pu',1.0, ...
		'angle', 0, ...
		'frequency', 60, ...
		'MVAsc3', 2000, ... 3?phase short circuit MVA= kVBase^2 / ZSC
		'MVAsc1', 2100, ...1?phase short circuit MVA. There is some ambiguity concerning the meaning of this quantity For the DSS, it is defined as kVBase^2 / Z1?phase where Z1?phase = 1/3 (2Z1+Z0) Thus, unless a neutral reactor is used, it should be a number on the same order of magnitude as Mvasc3.
		'x1r1', 4.0, ...
		'x0r0', 3.0, ...
		'Isc3', 10000, ...
		'Isc1', 10500, ...
		'R1', 1.65, ...
		'X1', 6.6, ...
		'R0', 1.9, ...
		'X0', 5.7, ...
		'ScanType','pos', ...
		'Sequence','pos', ...
		'Spectrum','defaultvsource', ...
		'BaseFreq', 60);
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