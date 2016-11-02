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

% get field names that need special handling (those that might contain
% several instances)
if ~isempty(s.dependency)
	dfn = fieldnames(s.dependency);
	% get all dependent vars (both parent and child)
	fn = dfn;
	for j = 1:length(dfn)
		fn = unique([fn s.dependency.(dfn{j})]);
	end
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
		case 'nphases'
			val{i} = dataclean(val{i},'phase');
		% Emergamps defaults to 1.35x normamps
		case {'name','wire'}
			val{i} = dataclean(val{i},'name');
		case 'like'
			if(ischar(val{i}))
				val{i} = dataclean(val{i},'name');
			elseif(isa(val{i},mfilename('class')))
				s.defaults = get(val{i});
				val{i} = s.defaults.Name;
			end
		case 'normamps'
			val{i} = dataclean(val{i},'num');
			s.defaults.Emergamps = val{i}*1.35;
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
	if ismember(propn{i},fn)
		% find out if prop is parent or child 
		if ismember(propn{i},dfn)
			% parent: increase parent index
			if length(s.data.(propn{i}))==1 && isempty(s.data.(propn{i}){1})
				id = 1;
			else
				id = length(s.data.(propn{i})) + 1;
			end
		else
			% child: find parent
			for k = 1:length(dfn)
				if ismember(propn{i},s.dependency.(dfn{k}))
					break;
				end
			end
			% child: get parent index
			id = length(s.data.(dfn{k}));
		end

		% set prop's value to appropriate index
		if id > 0
			s.data.(propn{i}){id} = val{i};
		else
			s.data.(propn{i}){1} = val{i};
		end
	else
		s.data.(propn{i}) = val{i};
	end
end

end
