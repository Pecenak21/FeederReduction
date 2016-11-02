function obj = dsstshape(s,varargin)
% OpenDSS TShape object
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
		'npts',0,... %Max number of points to expect in temperature shape vectors. This gets reset to the number of Temperature values found if less than specified.
		'interval',1,...%Switch input of  temperature curve data to a csv file containing (hour, Temp) points, or simply (Temp) values for fixed time interval data, one per line. NOTE: This action may reset the number of points to a lower value.
		'temp',[],...%Array of temperature values.  Units should be compatible with the object using the data. You can also use the syntax: Temp = (file=filename)     !for text file one value per line;Temp = (dblfile=filename)  !for packed file of doubles;Temp = (sngfile=filename)  !for packed file of singles. Note: this property will reset Npts if the  number of values in the files are fewer.
		'hour',[],...%Array of hour values. Only necessary to define this property for variable interval data. If the data are fixed interval, do not use this property. You can also use the syntax: hour = (file=filename)     !for text file one value per line; hour = (dblfile=filename)  !for packed file of doubles; hour = (sngfile=filename)  !for packed file of singles 
		'mean',0,...%Mean of the temperature curve values.  This is computed on demand the first time a value is needed.  However, you may set it to another value independently. Used for Monte Carlo load simulations.
		'stddev',0,...%Standard deviation of the temperatures.  This is computed on demand the first time a value is needed.  However, you may set it to another value independently.Is overwritten if you subsequently read in a curve. Used for Monte Carlo load simulations.
		'csvfile','',...%Switch input of  temperature curve data to a csv file containing (hour, Temp) points, or simply (Temp) values for fixed time interval data, one per line. NOTE: This action may reset the number of points to a lower value.
		'sngfile','',...%Switch input of  temperature curve data to a binary file of singles containing (hour, Temp) points, or simply (Temp) values for fixed time interval data, packed one after another. NOTE: This action may reset the number of points to a lower value.
		'dblfile','',...%Switch input of  temperature curve data to a binary file of doubles containing (hour, Temp) points, or simply (Temp) values for fixed time interval data, packed one after another. NOTE: This action may reset the number of points to a lower value.
		'sinterval',3600,...%Specify fixed interval in SECONDS. Alternate way to specify Interval property.
		'minterval',60,...%Specify fixed interval in MINUTES. Alternate way to specify Interval property.
		'action','',...%{DblSave | SngSave} After defining temperature curve data... Setting action=DblSave or SngSave will cause the present "Temp" values to be written to either a packed file of double or single. The filename is the Tshape name. 
		'Like','');
	obj.fieldnames = fieldnames(obj.defaults);
	% set the data structure to be a blank version of the defaults
	% structure
	obj.data = obj.defaults;
	obj.data(1) = []; obj.data(1).Name = ''; obj.data(1).Like = '';
	% some fields require re-naming
	obj.namemap = struct('R','%R','Rs','%Rs');
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