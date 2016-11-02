function [out, d, l, r] = structcmp( a , b )
%STRUCTCMP quickly compare two structs
%   STURCTCMP takes two structs as arguments and prints out any fields on
%   which the two differ

if ~isstruct(a) || ~isstruct(b)
	error('structcmp:invalidinput','Input is not struct');
end
if(nargout >=1 )
	d = {};
	l = {};
	r = {};
end

afn = fieldnames(a);
bfn = fieldnames(b);

if length(afn) ~= length(bfn)
	if(nargout == 0)
		fprintf('\n2 structs have different number of fields.');
	end
	id = ~ismember(afn,bfn);
	if sum(id) > 0
		if(nargout>=2)
			l = afn(id);
		else
			fprintf('\n  Struct A has following fields that B doesn''t:');
			for i = find(id)
				fprintf('  %s',afn{i});
			end
		end
	end
	
	id = ~ismember(bfn,afn);
	if sum(id) > 0
		if(nargout>=3)
			r = bfn(id);
		else
			fprintf('\n  Struct B has following fields that A doesn''t:');
			for i = find(id)
				fprintf('  %s',bfn{i});
			end
		end
	end
end

if(nargout == 0)
	fprintf('\n\nOther diffences:\n');
end
for i=1:length(afn);
	if ~isequaln(a.(afn{i}),b.(afn{i}))
		if(nargout >= 1)
			d = [d afn{i}]; %#ok<AGROW>
		else
			fprintf('''%s'' differs\n',afn{i});
		end
	end
end

if isempty(d) && isempty(l) && isempty(r)
    out = 1;
else
    out = 0;
end
end

