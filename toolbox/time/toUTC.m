function t = toUTC(time,location,DSTflag)
% convert time to UTC based on location of interest and whether it uses daylight saving time (DST) 
%
% Example of use:
%				t = toUTC
%				t = toUTC(now,'Los Angeles',0)

if nargin < 1
	time = now;
end
if nargin < 2
	location = 'Los Angeles';
end
if nargin < 3
	DSTflag = 0;
end

persistent tz;
if isempty(tz)
	tz = timeZones();
end

if DSTflag 
	t = tz.dst2utc(time,location);
else
	t = tz.st2utc(time,location);
end

end