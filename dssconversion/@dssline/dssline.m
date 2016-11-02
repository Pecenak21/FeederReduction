function obj = dssline(s,varargin)
% OpenDSS line object
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
	obj.defaults = struct( ...
		'Name','',...
		'Geometry','',...
		'LineCode','',...
		'Phases',3,...
		'bus1','',...
		'bus2','',...
		'Units','None',...	% Units can be used to allow specification of impedance and line length in different units.  See below.
		... % we want to make sure that units are printed before all of the quantities they apply to, so that the values won't accidentally get converted to a new system of units when we read in the file
		'Length',1,... %kft
		'Switch','no',... % switch effects all the R properties, so we place it before them so that if they are set too, they will be written after and not get lost
		'R1',0.058,... %ohms per 1000 ft
		'X1',0.1206,... %ohms per 1000 ft
		'R0',0.1784,... %ohms per 1000 ft
		'X0',0.4047,... %ohms per 1000 ft
		'C1',3.4,... %nF per 1000 ft
		'C0',1.6,... %nF per 1000 ft
		'BaseFreq'  ,60     ,...
		'Normamps'  ,[]     ,...
		'Emergamps' ,[]     ,...
		'Faultrate' ,0.0005     ,...
		'Pctperm'   ,[]     ,...
		'Rg'        ,0.01805     ,...
		'Xg'        ,0.155081     ,...
		'Rho'       ,100     ,...
		'Like'      ,''     ,...
		'Repair',[],... %hours to repair
		'Rmatrix',[],...
		'Xmatrix',[],...
		'Cmatrix',[],...
		'EarthModel','',...
		'enabled','' ...
		);
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