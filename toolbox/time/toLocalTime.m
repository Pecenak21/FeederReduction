function t = toLocalTime(UTCtime,location,DSTflag)
% convert UTC time to local time based on location of interest and whether it uses daylight saving time (DST) 
%
% Example of use:
%				t = toLocalTime
%				t = toLocalTime(nowUTC,'Los Angeles',0)

if nargin < 1
	UTCtime = nowUTC;
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
	t = tz.utc2dst(UTCtime,location);
else
	t = tz.utc2st(UTCtime,location);
end

end