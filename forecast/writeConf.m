function ok = writeConf(s, fn, overwrite)
% writeConf writes conf files that are sometimes readable by readConf
%
% It's not been thoroughly tested, but I needed a way to output some
% structs quickly, and matlab didn't want to write them for me, so I threw
% this together quickly.  Contributions are welcome.
%
% Usage: writeconf( structdata, filename, [overwrite_flag])
%
% overwrite_flag = true must be specified to overwrite an existing file.
%
% See Also: readConf

ok = false;

% check inputs (file exist? struct data?)
if(exist(fn,'file') && (~exist('overwrite','var') || ~overwrite))
	warning('writeConf:fileExists','defaulting to not overwriting an existing conf file');
	return;
end
if(~isstruct(s))
	error('writeConf:needStruct','writeConf only works with structs');
end

% open file
fid = fopen(fn,'w');
if(fid==-1), return; end

% print the fields, one per line
fn = fieldnames(s);
for i=1:numel(fn)
	if(isnumeric(s.(fn{i})))
		fprintf(fid,'%s\t%g\n',fn{i}, s.(fn{i}));
	elseif iscell(s.(fn{i}))
		% cells print as
		% $$FIELDNAME
		% Value1
		% Value2
		% ...
		% $$
		fprintf(fid,'$$%s\n',fn{i});
		x = s.(fn{i});
		for j=1:numel(x)
			if(isnumeric(x{j}))
				fprintf(fid,'%g\n',x{j});
			else
				fprintf(fid,'%s\n',x{j});
			end
		end
		fprintf(fid,'$$\n');
	else
		fprintf(fid,'%s\t%s\n',fn{i}, s.(fn{i}));
	end
end

% close file, cleanup output
ok = ~fclose(fid);
if(~nargout)
	clear('ok');
end

end