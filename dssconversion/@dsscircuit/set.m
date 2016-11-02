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
	% as we go through the main body of the loop, we'll tweak the value
	% itself, and make any necessary changes to OTHER properties, and then
	% let one line be responsible for setting the current property
	switch(lower(propn{i}))
		case {'name','bus2','spectrum','like'}
            val{i} = dataclean(val{i},'name');
		case 'bus1'
			val{i} = dataclean(val{i},'name');
			s.defaults.bus2 = [val{i} '.0.0.0'];
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