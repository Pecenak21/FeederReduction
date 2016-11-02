function [sr ss] = sunriseSunset(position, day)
% sunriseSunset(position, day)
% look up the sunrise and sunset times for the day in question
%	sunrise and sunset times are found by looking for the solar zenith
%	angle to be less than 93 degrees, and are rounded to the half-minute
%
% This function is designed to determine the first and last timestamps that
% will be captured for a given imager for the day

% get day in matlab datenum format for the start of the day (in UTC)
if(ischar(day) || iscell(day))
	day = datenum(day);
end
tz = round(position.longitude/15); % 15 degrees per hour
day = floor(day);

% Attempt to calculate using the c-based API with direct lookup of sr/ss
% if that fails, we fallback to running siSunPosition for the whole day.

% The fast method is currently disabled because it lacks the flexibility to get the exact times that the USI starts/stops capturing (93 degree SZA), and since we're not using this in time-dependent code anyway.

try
	[~,~,~,sr,ss] = siSunPosition_mex(day+0.5, position, tz);
	% the C solar position algorithm returns fractional hours instead of an absolute time.
	ss = (ss-tz)/24 + day(:);
	sr = (sr-tz)/24 + day(:);
	% In some locations, the SPA seems to return the sunset time of the
	% previous day instead of the sunset time of the current day, but it's
	% usually still to within the 30 seconds, so I'm not going to worry
	% about it for now
	
	%return;
catch e
	% we ignore errors that would result from the mex code not having been compiled
	if(~strcmp(e.identifier,'MATLAB:UndefinedFunction'))
		rethrow(e);
	end
end

% adjust to have approximately local midnight:
day = day - tz/24;

% allocate outputs
sr = zeros(size(day));
ss = sr;

% calculate the solar angles all day each day, and then compare to get 
for i = 1:length(day);
	day_ = day(i);
	tr = day_:1/24/2:day_+1;
	for j = 1:2 % increase precision an extra time
		p = siSunPosition(tr, position);
		isDaytime = [p.zenith]<93;
		si = find(isDaytime, 1, 'first');
		ei = find(isDaytime, 1, 'last');
		sv = tr(si); ev = tr(ei);
		dt = tr(2)-tr(1);
		tr = [tr(si-1):dt/60:sv, ev:dt/60:tr(ei+1)];
	end
	sr(i) = sv;
	ss(i) = ev;
end

end
