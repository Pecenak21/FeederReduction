function v = get(s,propn)
% the get function for dsscapacitor objects is responsible for filling in
% default values so that a user or other code will be aware of what OpenDSS
% does

% start by copying the structure
v = reshape([s.data], size(s));

% patch up the structure to fill in defaults
if(nargin > 1)
	fn = propn;
	if(~iscell(fn)), fn = {fn}; end
elseif(length(v)>=1)
	fn = s(1).fieldnames;
else
	fn = {};
end
for i=1:length(fn)
	for j=1:numel(s)
		% match field name (case-insensitive)
		try 
			v(j).(fn{i});
		catch
			fns = fieldnames(s);
			[id,id] = ismember(lower(fn{i}),lower(fns));
			fn{i} = fns{id}; 
			if ~iscell(propn)
				propn = fns{id};
			else
				propn{i} = fns{id};
			end
		end
		if(isempty(v(j).(fn{i})))
			v(j).(fn{i}) = s(j).defaults.(fn{i});
		end
	end
end

% return only a part of it if requested
if(nargin > 1 && ~isempty(v))
	if(~iscell(propn)), propn = {propn}; end
	for i=1:numel(propn)
		v_(:,i) = reshape({v.(propn{i})},[numel(s),1]);
	end
	v = v_;
	if(numel(v) == 1)
		v = v{1};
	end
end
end
