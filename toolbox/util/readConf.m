function [ conf ] = readConf( conf_path , convertFlag )
%READCONF Reads a configuration file
%   readConf is intended as a native alternative to
%   bu.util.conf.Configuration, with the intent of easier modifications and
%   eventual extension to updating conf files as well.
%
% Usage:
%	conf = readConf( conf_path )
%	conf will be a matlab structure of strings
%	conf_path should be the path to the file

conf_path = char(conf_path); % handle java.lang.String and java.io.file inputs

fid = fopen(conf_path);
if (fid == -1)
	error('Could not open conf file %s\nPerhaps it doesn''t exist?',conf_path);
end

if ~exist('convertFlag','var')
	convertFlag = 1;
end

% return an empty struct even for empty conf files
conf = struct();
% readMulti is a flag for determining whether we're in the middle of a
% multi-line sequence right now
readMulti = false;

tline = fgetl(fid);
while ischar(tline)
	% skip comment lines and ignore inline comment
    id = find(tline=='#',1);
    if(~isempty(id)), tline(id:end) = []; end;
	% trimming space before and after a line
	tline = strtrim(tline);
    if isempty(tline); tline = fgetl(fid); continue; end;
	
	% special handling for multivalue key/value set
	if( length(tline)> 1 && all(tline(1:2)=='$') )
		% try to read the key
		tok = regexp(tline,'^\$\$\s*(\w+)\s*$','tokens');
		if ~isempty(tok)
			tok = tok{1};
			multiKey = tok{1};
			try
				% check if key is already in use
				if isfield(conf,multiKey)
					warning('readConf:multi','Multiple copies of key %s; will be using the last value specified',multiKey);
				end
				% create an empty cell array to save multiple values in
				conf.(multiKey) = {};
			catch err
				warning('readConf:badField','Config file contains a parameter that is invalid as a structure field:\n%s',tok{1});
			end
			readMulti = true;
			multiN = 1;
		else
			readMulti = false;
		end
	elseif readMulti % in a multi-line read
		if ~any(tline(1:min(length(tline),2))=='%') %% lines are worthless; I'm ignoring them, but I think they might be important to Bryan's Java code
			% read the value
			tok = regexp(tline,'^\s*(.+)\s*$','tokens');
			if(~isempty(tok))
				tok = tok{1};
				conf.(multiKey){multiN} = tok{1};
				multiN = multiN+1;
			end
		end
	elseif( length(tline)>1 && all(tline(1:2)=='%') )
		warning('readConf:corrupt','Found %% outside multiline value; either the file is corrupt or you accidentally tried to use matlab-style comments in the conf file!');
	else % regular line
		% split the line into a key and a value
		tok = regexp(tline,'^(\w+)\s+(.+)\s*$','tokens');
		if(~isempty(tok))
			tok = tok{1};
			try
				% check if key is already in use
				if isfield(conf,tok{1})
					warning('readConf:multi','Multiple copies of key %s; will be using the last value specified',tok{1});
				end
				% save value
				conf.(tok{1}) = tok{2};
			catch err
				warning('readConf:badField','Config file contains a parameter that is invalid as a structure field:\n%s',tok{1});
			end
		end
	end
	tline = fgetl(fid);
end

fclose(fid);

% Convert all number fields with strings in conf to numbers
if convertFlag
	conf = parseConf(conf);
end

end

