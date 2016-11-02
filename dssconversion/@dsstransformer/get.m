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
% handle filling in extra data around dependencies
try
	s.dependency;
	hasdep = true;
catch
	hasdep = false;
end
if(hasdep && numel(v)>=1)
	fns = find(ismember(s(1).dependency(2:end),fn));
	fns = s(1).dependency(fns+1);
	for i=1:length(fns)
		for j=1:numel(s)
			% if the data array is too short (this happens when we haven't
			% explicitly set a number of windings), extend it
			m = length(s(j).defaults.(fns{i}));
			if(length(v(j).(fns{i})) < m)
				v(j).(fns{i}){m} = [];
			end
			% fill in missing elements
			m = cellfun(@isempty,v(j).(fns{i}));
			v(j).(fns{i})(m) = s(j).defaults.(fns{i})(m);
		end
	end
end

% Handle plural properties
% This block only effects the transformer class
if(strcmp(mfilename('class'),'dsstransformer'))
	fns = {'kVs','kVAs','Taps','Rs','Buses','Conns';...
		'kV', 'kVA', 'tap', 'R','Bus','Conn'};
	fns_ = ismember(fns(1,:),fn);
	for i=find(fns_)
		for j=1:numel(s)
			l = length(v(j).(fns{1,i}));
			if(i>=5) %i.e. for fns {buses, conns}
				v(j).(fns{1,i}) = [v(j).(fns{1,i}) s(j).defaults.(fns{2,i})(l+1:end)];
			else
				v(j).(fns{1,i}) = [v(j).(fns{1,i}) s(j).defaults.(fns{2,i}){l+1:end}];
			end
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
