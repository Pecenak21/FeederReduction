function obj = dssgenerator(s,varargin)
% OpenDSS generator object
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
	% PLEASE refer to OpenDSS manual Program Revision: 7.4.3 page 125 for
	% more details about each properties
	obj.defaults = struct( ...
        'Name','',...
		'Phases',3,...
		'bus1','',...% Name of bus to which the generator is connected. Include node definitions if the terminal conductors are connected unusually. 3?phase Wye?connected generators have 4 conductors; Delta?connected have 3. Wye?connected generators, in general, have one more conductor than phases. 1?phase Delta has 2 conductors; 2? phase has 3. The remaining Delta, or line?line, connections have the same number of conductors as phases.
		'Kv',12.47,... %Base voltage for generator. For 2? or 3?phase generators, specified in phase?tophase kV. For all other generators, the actual kV across the generator branch. If wye (star) connected, specify the phase?to?neutral (L?N) kV. If delta or phase?tophase connected, specify the phase?to?phase (L?L) kV.
        'Kw',[],...% nominal kW for generator. Total of all phases.
        'Pf',[],...
        'Model',1,...
        'Yearly',{[]},...
        'Daily',{[]},...
        'Duty',{[]},...
		'Dispvalue',0.0,...
        'Conn','wye',...
        'Kvar',[],...
        'Rneut',0,...
        'Xneut',[],...
        'Status','variable',...
        'Class',1,...
		'Maxkvar',[],...
		'Minkvar',[],...
		'Pvfactor',0.1,...
		'Debugtrace','No',...
		'ForceON','No',...
        'Vminpu',0.95,...
        'Vmaxpu',1.05,...
        'kVA',[],...
		'MVA',[],...
		'Xd',1.0,...
		'Xdp',0.27,...
		'Xdpp',0.2,...
		'H',1.0,...
		'D',1.0,...
		'UserModel','',...
		'UserData','',...
		'ShaftModel','',...
		'ShaftData','',...
		'spectrum','default',...
        'BaseFreq',60,...%Base frequency for which this generator is defined. Default is 60.0.
        'Like'      ,'',...
		'Enabled','true' ...%{yes|no or true|false}
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