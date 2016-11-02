function s = fieldnamefix(s,validnames)
% fieldnamefix(sin,validnames) rewrites the field names in the input
%	structure sin to be on the list 'validnames'.  Fields are duplicated
%	rather than removed.

%check through each field name
waschar = false;
if(isstruct(s))
	fn = fieldnames(s);
elseif(iscellstr(s))
	fn = s;
elseif(ischar(s))
	fn = {s};
	s = {s};
	waschar = true;
else
	error('first input must be a struct or a cell array of strings');
end
for i=1:length(fn);
	% if the field name is already valid, do nothing
	match = strcmpi(fn{i},validnames);
	if(any(match))
		if(any(strcmp(fn{i},validnames(match)))), continue; end
	else
		% check for matches up to the length of the fieldname we're testing
		match = strncmpi(fn{i},validnames,length(fn{i}));
		% if no matches, give up
		if(~any(match))
			if(~isstruct(s)), s{i} = ''; end
			continue;
		end
	end
	newn = validnames{find(match,1)};
	% don't overwrite existing correct fields
	if(isstruct(s))
		if(isfield(s,'newn')), continue; end
		s.(newn) = s.(fn{i});
	else
		s{i} = newn;
	end
end
if(waschar)
	s = s{1};
end

end
