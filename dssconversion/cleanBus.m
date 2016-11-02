function o = cleanBus(x)
% clean bus id from 1231231.1.2 to 1231231

if iscell(x)
	for i = 1:length(x)
		o{i} = stripPhase(x{i});
	end
else
	o = stripPhase(x);
end

end

function o = stripPhase(x)
	% find the dot
	id = find(x=='.'); 
	if ~isempty(id)
		o = x(1:id(1)-1);
	else
		o = x;
	end
end

