function obj = dsspvsystem(s,varargin)
% OpenDSS PVsystem object
%
% To implement PVsystem using forecasted/ measured power output,  use 'pmpp' 
% as rated PV output (in kW), irradiance=1, PF=1, kVA (inverter rating) > pmpp, 
% and use loadshape to define forecasted/real daily or yearly PV output 
% with values between 0 and 1.
%
% To implement PVSystem object properly using irrad and temp info, 3 important basic parameters are needed:
%    1. Average "Pmpp" for panel at 1kW/m^2 irradiance and 25C.
%    2. Per unit variation "PTcurve" of Pmpp vs Temp at 1kW/m^2 irradiance.
%    3. Efficiency curve for the inverter "EffCurve", per unit efficiency vs per unit power.
% The panel output is estimated by: 
%    Panel kW = Pmpp (in kW @1kW/m2 and 25 C) * Irradiance (in kW/m2) * Factor(@actual T)
% For more details on how to model PV system, please read  <a href="matlab:web('http://svn.code.sf.net/p/electricdss/code/trunk/Distrib/Doc/OpenDSS%20PVSystem%20Model.pdf')">this document</a>.
%
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
		'phases',3,...
        'bus1','',...
        'kv',12.47,...
        'irradiance',1,...
        'Pmpp',500,...
        'Temperature',25,...
        'pf',1,...
        'conn','wye',...
        'kvar',[],...
        'kVA',500,...
        'cutin',20,...
        'cutout',20,...
        'EffCurve','',...
        'PTCurve','',...
        'R',0,...
        'X',50,...
        'model',1,...
        'Vminpu',.9,...
        'Vmaxpu',1.1,...
        'yearly','',...
		'daily','',...
        'duty','',...
        'Tyearly','',...
        'Tdaily','',...
        'Tduty','',...
        'class',1,...
        'UserModel','',...
        'UserData','',...
        'debugtrace','No',...
        'spectrum','default',...
        'basefreq',60,...
        'enabled','Yes',...
        'Like',{''}...
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
    % some fields require re-naming
	obj.namemap = struct('cutin','%cutin','cutout','%cutout','R','%R','X','%X','PTCurve','P-TCurve');
end

% make it an object
obj = class(obj,cn);

% special messgage when object is created
persistent warningOn;
if isempty(warningOn), warningOn = 1;
else warningOn = 0;
end
% if warningOn
% 	warning(sprintf(['To implement PVSystem object properly. 3 important basic parameters are needed:\n'...
% 		'\t1. Average "Pmpp" for panel at 1kW/m^2 irradiance and 25C.\n'...
% 		'\t2. Per unit variation "PTcurve" of Pmpp vs Temp at 1kW/m^2 irradiance.\n'...
% 		'\t3. Efficiency curve for the inverter "EffCurve", per unit efficiency vs per unit power.\n'...
% 		'The panel output is estimated by: \n'...
% 		'\tPanel kW = Pmpp (in kW @1kW/m2 and 25 C) * Irradiance (in kW/m2) * Factor(@actual T)\n'...
% 		'\nFor more details on how to model PV system, please read  <a href="./doc/OpenDSS PVSystem Model.pdf">this document</a>.']));
% end

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