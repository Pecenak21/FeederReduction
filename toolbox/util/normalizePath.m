function p = normalizePath(filepath)
%normalizePath dereferences drive ids and converts relative paths to be absolute paths
% paths that begin with a '$' are expanded to dereference the following driveID.  As a special case, paths that begin with '$$' are dereferenced to point to the root of the source repository
%
% Drive IDs allow us to have a general pointer to a mounted drive across multiple systems; see getDriveById()
%
% For example:
% normalizePath('$KleisslLab4TB1') -> /mnt/lab_4tb1
% normalizePath('$$/local')        -> /path_to_code/conf/local
% normalizePath('local')           -> /path_to_code/conf/local  (after running def_addpath)
%
% default is to return the passed path

% accept cell input
if(iscell(filepath))
	p = cellfun(@normalizePath,filepath,'UniformOutput',0);
	return;
end

% default behavior:
p = filepath;

% if the first char is '$' expand a tag
if(filepath(1)=='$')
	% if the second char is also '$', replace with the source root
	if(filepath(2) == '$')
		filepath(1:2) = [];
		if(~isempty(filepath) && all(filepath(1)~='/\')), filepath = ['/' filepath]; end %take care not to insert a double slash
		p = [fileparts(which('def_addpath.m')) filepath];
	else % dereference the driveid
		en = min([length(filepath), find(filepath=='/',1)-1, find(filepath=='\',1)-1]);
		p = getDriveById(filepath(2:en));
		if(isempty(p))
			warning('siNormalizePath:NoDrive','Couldn''t find mount point with requested drive_id, %s',filepath(2:en));
			p = filepath; return;
		end
		filepath(1:en) = [];
		p = [p filepath];
	end
elseif isrelpath(filepath)
	if(~exist(filepath,'file'))
		warning('siNormalizePath:NoFile','''%s'' doesn''t exist! I can''t normalize a path I can''t find.',filepath);
		p = filepath;
		return;
	end
	if(~isempty(which(filepath))) % file is in the search path
		p = which(filepath);
	else
		try
			p = what(filepath);
			p = p(1).path;
		catch e %#ok<NASGU>
			p = java.io.File( pwd , filepath );
			p = char(p.getCanonicalPath());
		end
	end
end

%% Checks for relative paths
%  returns true if the path is relative, and false if it is not
	function f = isrelpath(ip)
		f = true;
		% unix paths (and maybe some absolute paths on windows?) start with /
		if(any(ip(1)=='/\'))
			f = false; return;
		end
		% windows paths can start with a letter (which I'm lazy about checking) followed by a ':'
		% and a slash (if we're not just referencing the root)
		if(ispc() && length(ip)>=2 && ip(2)==':' && (length(ip)<3 || any(ip(3)=='/\')))
			f = false; return;
		end
		
	end

if( ispc() )
	p = java.io.File(p);
	p = char(p.getCanonicalPath());
	
	% replace the stupid backslash in windows by forward slash
	p(p=='\') = '/';
end

end