function s = char(o)
% Convert object to string that can be put in an openDSS file
%	we are careful only to write values that were set by the user, not
%	defaults

%% handle array inputs by recursion
if(length(o)~=1) % actually handles the 0 case nicely too!
	s = '';
	for i=1:numel(o)
		s = [s char(o(i))];
	end
	return;
end

%% setup to export data fields
% we only export data, not defaults, and we'll do the Name field manually,
% so we start with the second fieldname.
dcap = o.data;
fn = o.fieldnames(2:end);
% skip any dependencies on the first pass
try
	fn = setdiff(fn,o.dependency);
catch
end

%% output 'name' field
% we also want to know the opendss class name, which is by convention the
% matlab class name with 'dss' left off (or in an earlier version 'D' left
% off
cname = mfilename('class');
if(strmatch('dss',cname)), cname = cname(4:end);
elseif(strmatch('D',cname)), cname = cname(2:end);
end
cname(1) = upper(cname(1));

% Add the name
s = ['New ' cname '.' dcap.Name];

%% Loop through and output any other properties as needed
try
	dorename = (~isempty(o.namemap));
catch
	dorename = 0;
end
for i=1:length(fn)
	% rename the field if appropriate
	if(dorename && isfield(o.namemap,fn{i}))
		fname = o.namemap.(fn{i});
	else
		fname = fn{i};
	end
	if ~isempty(dcap.(fn{i}))
		% grab the value in a variable to make the remaining lines a little
		% shorter
		s = outputstring(s,fname,dcap.(fn{i}));
	end
end
%now handle the dependencies if present
try
	fn = o.dependency;
catch
	fn = {};
end
if(~isempty(fn))
	% grab all the data values to make it easy to see which indices we need
	% to write
	dgrid = cell(1);
	for i=2:length(fn)
		d = dcap.(fn{i});
		if(isempty(d)), continue; end;
		dgrid(i-1,1:length(d)) = d;
	end
	% handle output name mapping
	for i=1:length(fn)
		try
			fn{i} = o.namemap.(fn{i});
		catch
		end
	end
	% iterate over indices that have data
	m = ~cellfun('isempty',dgrid);
	for i=find(any(m))
		s = sprintf('%s %s=%i',s,fn{1},i);
		for j=find(m(:,i))'
			s = outputstring(s,fn{j+1},dgrid{j,i});
		end
	end
end

% when we're all done, add a trailing newline so we're all ready to go into
% a file
s = sprintf('%s\n',s);

end

function s = outputstring(s,fn,v)
% Handle different cases of data type
if isempty(v)
	return;
end

% string data:
if(ischar(v))
	s = [s ' ' fn '=' v ''];
elseif(islogical(v)||isnumeric(v))
	% matrix or vector data
	% note that we don't have a good way to automatically detect
	% vectors that might have been intended to be matrices, so it's
	% important for the parsing code to interpret things like
	% [1 | 1 2 | 1 2 3] correctly.
	if(length(v)>1)
		if(length(v) ~= numel(v))
			% this appends a '| character at the end of each line of
			% the matrix and then wraps the string back onto one line.
			% e.g. [1 2 3| 4 5 6| 7 8 9]
			v = [num2str(v), repmat('| ',[size(v,1),1])];
			v = reshape(v',[1,numel(v)]);
			s = [s ' ' fn '=[' v(1:end-1) ']'];
		else
			s = [s ' ' fn '=[' num2str(v) ']'];
		end
	else %and finally scalar numeric data:
		s = [s ' ' fn '=' num2str(v)];
	end
elseif(iscell(v))
	% if v is an array of strings (for Buses or Conns)
	if length(v) == 1
		s = [s ' ' fn '=' v{1} ];
	else
		str = v{1};
		for j = 2:length(v)
			str = [str ' ' v{j}];
		end
		s = [s ' ' fn '=[' str ']'];
	end
else
	error('Do not know that kind of data!');
end
end