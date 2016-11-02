function obj = dssxycurve(s,varargin)
% OpenDSS XYCurve object
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
		'npts',0,... %Max number of points to expect in curve. This could get reset to the actual number of points defined if less than specified.
		'Points','',...% One way to enter the points in a curve. Enter x and y values as one array in the order [x1, y1, x2, y2, ...]. For example:
					...% Points=[1,100 2,200 3, 300] 
					...% Values separated by commas or white space. Zero fills arrays if insufficient number of values.
		'Xarray',[0],...%Alternate way to enter X values. Enter an array of X values corresponding to the Y values.  You can also use the syntax:  Xarray = (file=filename)     !for text file one value per line. Note: this property will reset Npts to a smaller value if the  number of values in the files are fewer.
		'Yarray',[0],...%Alternate way to enter Y values. 
		'x',0,... %Enter an value and then retrieve the interpolated Y value from the Y property.
		'y',0,... %Enter an value and then retrieve the interpolated X value from the X property.
		'csvfile','',...%Switch input of  X-Y curve data to a CSV file containing X, Y points one per line. NOTE: This action may reset the number of points to a lower value.
		'sngfile','',...%Switch input of  X-Y curve data to a binary file of SINGLES containing X, Y points packed one after another. NOTE: This action may reset the number of points to a lower value.
		'dblfile','',...%Switch input of  X-Y  curve data to a binary file of DOUBLES containing X, Y points packed one after another. NOTE: This action may reset the number of points to a lower value.
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