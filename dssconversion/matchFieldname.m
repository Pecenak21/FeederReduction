function [ fname ] = matchFieldname( fname, target )
% matching a field name to its target structure (making things case
% insensitive)
% input:
%			target: is either a struct or a cell of field names. If a struct then the fieldnames are extracted automatically 

if istruct(target) 
	fns = fieldnames(target);
elseif iscell(target)
	fns = target;
else
	error('Not supported target type. Must be a cell of strings or a struct!');
end

[~, id] = ismember(lower(fname),fns);
if ~isempty(id)
	fname = fns{id};
else
	warning('Target doesn''t include specified field! Check the fieldname again!');
end

end

