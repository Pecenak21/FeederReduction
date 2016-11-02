function sun = siSunPosition(time, latitude, longitude, altitude)
% siSunPosition(time, latitude, longitude, altitude) looks up the solar position using NREL's Solar Position Algorithm (SPA).
% siSunPosition(time, position)
%
% Inputs:
%	* time - one or more matlab datenums describing the times for which to look up the solar position
%	* latitude in degrees N
%	* longitude in degrees E
%	* altitude of the site in meters
%
%	Alternatively, pass a struct with fields for latitude, longitude, and altitude as the second argument
%
% Output:
%	An array of structs with the following fields: time, zenith, azimuth, earthsundistance

persistent doCompile;

if( isstruct(latitude) || isa(latitude,'bu.science.geography.Position') )
	pos = latitude;
else
	pos = struct('latitude',latitude,'longitude',longitude,'altitude',altitude);
end

%% First try to compile and use the C version of this routine
%  Originally I found a nifty 'autocompile' function that you could drop into the mfile and basically it would compile the c code into a .mex file if needed, and run that, and after that, the mex version is always used instead.
%  The problem with that is that really we only care about the C as an interface to the existing spa.c (which it is desirable to be able to use directly from NREL; otherwise we'd just port that to matlab and be done)
%  But there are some operations (for example, checking if the compiled .mex is out of date, or fetching a table of leap seconds) that it might be nice to be able to do, and that we currently don't do.  These will be _much_ easier to write in matlab than in C, so it's nice to have a wrapper .m that's always called so that we can handle those if we ever feel the need.
if( (isempty(doCompile) && needsCompile('siSunPosition_mex') ) || doCompile)
	compileHelper('siSunPosition_mex.c',{'spa.c'});
	doCompile = 0;
end
if(~isempty(which('siSunPosition_mex')))
	% call the mex helper function
	[zenith, azimuth, earthsundistance] = siSunPosition_mex(time,pos);
	sun = struct('time',time,'zenith',zenith,'azimuth',azimuth,'earthsundistance',earthsundistance);
	return;
end

%% If the plan to use the C version falls through, try java:
sun = siSolarAngles(bu.util.Time.datevecToTime(datevec(time)), pos);

% replace buTime with matlab time
for i=1:length(sun)
	sun.time = time(:);
end

	%% needsCompile
	%  needsCompile does a quick check to see whether the mex file exists, and another to see if it's newer than the source file
	%  right now it's only actually called once, so the check for an existing compiled version is less useful than it might be, but this avoids the extra overhead of needing to check every time
	function flag = needsCompile(fn)
		mc = which(fn);
		if(isempty(mc))
			flag = true; return;
		end
		% this implies that we found the compiled version
		mi = which([fn '.c']);
		mi = dir(mi); mc = dir(mc);
		if(datenum(mi.date)>datenum(mc.date))
			doCompile = 1;
			flag = true;
		else
			doCompile = 0;
			flag = false;
		end
	end
end

function worked = compileHelper(src,deps)
%given paths are relative to the source file
mexdir = fileparts(mfilename('fullpath'));
deps = strcat(mexdir, '/', deps);
if(~iscell(deps)), deps = {deps}; end

% try to compile the mex file on the fly
warning('siMEX:compiling','trying to compile MEX file from %s', src);
% for some stupid reason, MATLAB defaults to ANSI C on unix platforms, which means a lot of the comments in the SPA code don't work
if(isunix()), deps{end+1} = 'CFLAGS="\$CFLAGS\ -std=c99"'; end
% here's the compile:
worked = ~mex([mexdir '/' src],deps{:});

% try to move it to the right folder.  ignore failures, since they don't _really_ matter.
if(worked)
	[~] = movefile(regexprep(src,'\.c$',['.' mexext]),mexdir);
end

end
