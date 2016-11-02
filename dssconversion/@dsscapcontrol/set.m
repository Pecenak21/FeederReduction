function s = set(s,propn,val,varargin)
% The Set function is responsible for handling all the magic of making sure
% we format the object data correctly.  For example:
%	* Some kinds of data we try to correct subtle mistakes
%	* Setting some properties in opendss triggers other properties to be
%	cleared or have a default value
% All of this kind of detail is taken care of here in the set function

% handle cell value
while iscell(val) && length(val)==1
	val = val{1};
end
while iscell(propn) && length(propn)==1
	propn = propn{1};
end

% handle setting multiple values at once
if(nargin > 3)
	val = reshape([{propn,val} varargin],[2 (nargin-1)/2])';
	propn = fieldnamefix(val(:,1),s.fieldnames);
	val = val(:,2);
else %otherwise just make sure we're setting a valid value
	propn = fieldnamefix({propn},s.fieldnames);
	val = {val};
end

% set the data, altering a few cases as desired and setting other values or
% defaults in some other cases
for i=1:length(propn)
	if(isempty(propn{i})), continue; end
	% as we go through the main body of the loop, we'll tweak the value
	% itself, and make any necessary changes to OTHER properties, and then
	% let one line be responsible for setting the current property
	switch(lower(propn{i}))
		case {'eventlog','voltoverride'} %logical type, output as string/char
			val{i} = dataclean(val{i},'logical','string');
		case {'name','like','element','capacitor'}
			val{i} = dataclean(val{i},'name');
		case 'type'
			val{i} = controlTypeClean(val{i});
		case {'ptphase','ctphase'}
			val{i} = dataclean(val{i},'monitoredPhase');
		otherwise
			% determine property's type/class and convert given value to that type
			cl = class(s.defaults.(propn{i}));
			if ~isa(val{i},cl)
				switch cl
					case 'char'
						val{i} = char(val{i});
					case 'double'
						val{i} = dataclean(val{i},'num');
				end
			end
	end
	s.data.(propn{i}) = val{i};
end

end

function type = controlTypeClean(s)
%Clean up control type
if ~ischar(s)
	error('Invalid control type')
end

s = regexprep(s,'[\s"'']','');
switch lower(s(1))
	case 'v'
		type = 'voltage';
	case {'c','i'}
		type = 'current';
	case 'k'
		type = 'kvar';
	case 't'
		type = 'time';
    case 'p'
        type = 'PF';
end

end

function p = mornitoredPhaseClean(s)
% clean up mornitored phase.
% output:
%			p : {1,2,3} phase number; for delta/ll connection use the first of the two phases (so 1 for 1-2, 2 for 2-3, 3 for 3-1). Default to 1 if input is invalid.

if ischar(s)
	s = sort(lower(s));
	allphase = 'abcxyz';
	[val id] = ismember(s,allphase);
	ps = allphase(id);
	ps(ps=='a') = 'x';
	ps(ps=='b') = 'y';
	ps(ps=='c') = 'z';
	ps = unique(ps);
	switch ps
		case {'x','xy'}
			p = 1;
		case {'y','yz'}
			p = 2;
		case {'z','xz'}
			p = 3;
		otherwise
			warning('Invalid monitored phase input');
			p = 1;
	end
elseif isnumeric(s)
	p = uint8(s);
else
	warning('Invalid monitored phase input');
	p = 1;
end

end