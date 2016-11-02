function [ c ] = cleanDSScircuit( c )
% clean up DSSCircuit by: 1) fix transformer name for regcontrol (referring to
% transformer's bus)

if isfield(c,'regcontrol')
	for i = 1:length(c.regcontrol)
		[y, idx] = ismember( lower(c.regcontrol(i).transformer), lower({c.transformer.Name}) );
		if ~isempty(idx)
			c.regcontrol(i).transformer = c.transformer(idx).Name;
		end
	end
end

end

