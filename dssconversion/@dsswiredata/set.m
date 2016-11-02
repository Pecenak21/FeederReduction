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
		% some programs specify phases individually, but we just want the number
		case 'rac'
			val{i} = dataclean(val{i},'num');
			if isempty(s.data.Rdc)
				s.defaults.Rdc = val{i};
			end
		case 'rdc'
			val{i} = dataclean(val{i},'num');
			if isempty(s.data.Rac)
				s.defaults.Rac = val{i};
			end
		case 'radius'
			val{i} = dataclean(val{i},'num');
			if isempty(s.data.GMRac)
				s.defaults.GMRac = val{i}*.7788;
			end
		case 'gmrac'
			val{i} = dataclean(val{i},'num');
			if isempty(s.data.Radius)
				s.defaults.Radius = val{i}/.7788;
			end
		case 'normamps'
			val{i} = dataclean(val{i},'num');
			if isempty(s.data.Emergamps)
				s.defaults.Emergamps = val{i};
			end
		case 'emergamps'
			val{i} = dataclean(val{i},'num');
			if isempty(s.data.Normamps)
				s.defaults.Normamps = val{i};
			end
		case 'name'
			val{i} = dataclean(val{i},'name');
		case 'like'
			if(ischar(val{i}))
				val{i} = dataclean(val{i},'name');
			elseif(isa(val{i},mfilename('class')))
				s.defaults = get(val{i});
				val{i} = s.defaults.Name;
			end
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
