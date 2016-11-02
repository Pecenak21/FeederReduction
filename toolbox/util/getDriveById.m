function p = getDriveById(id)
% This function is intended as a cross platform method of identifying network drives across multiple computers
%
% It looks for a file called 'drive_id' at the root of every mount point on the system and returns the mount point that matches the requested id
%
% IDs are _NOT_ case sensitive on Unix/Mac. On Windows, it simply calls bu.io.File_bu.getDriveLetter, which I don't know anything about the internals of
% Note that although the drive ID comparison is not case sensitive, your filesystem may be, so you should make sure the file name is 'drive_id' all lowercase
%
% This function is similar to bu.io.File_bu.getDriveLetter, but works on Mac OS as well

if( isunix() )
	% mount command lists all mounted filesystems
	[~,x] = system('mount');
	% each line on a mac is:
	% ... on /path/to/mount/point (fstype,opts)
	% or on ubuntu linux:
	% ... on /path/to/mount/point 
	st = strfind(x,' on ')+4;
	en = strfind(x,' type ')-1;
	if(length(st) ~= length(en))
		en = strfind(x,' (')-1;
	end
	if(length(st) ~= length(en))
		error('couldn''t understand output from the mount command');
	end
	% look for a drive_id file on each mount point
	for i=1:length(st)
		p = x(st(i):en(i));
		fp = [p '/drive_id'];
		% if we find the file, read it
		if(exist(fp))
			fp = fopen(fp);
			l = fgetl(fp);
			fclose(fp);
			% if the file's contents match the requested id, return that mount point
			if(strcmpi(id,l)), return; end;
		end
		% if we didn't get a match, set the path empty so that if we exit on this loop
		p = [];
	end
end
% otherwise, ispc()
p = char(bu.io.File_bu.getDriveLetter(id));
if(isempty(p))
	p = char(bu.io.File_bu.getDriveLetter(upper(id)));
end

end
