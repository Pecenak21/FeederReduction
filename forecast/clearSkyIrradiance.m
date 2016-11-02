function csk = clearSkyIrradiance( pos, time, tilt, azimuth )
% This function calculates the clear sky GI, DNI and DI from Linke Turbidity database
% If tilt and azimuth is not specified, they are assumed to be 0 (flat) and 180 (south) respectively
%
% Input:
%		pos		is a struct, specifying latitude, longitude, altitude (meters)
%				(bu.science.geography.Position should work too)
%		time	(datenum) at which to perform the calculation
%		tilt	(degrees) of the solar array; flat is 0
%		azimuth (degrees) of the solar array; N is 0, E is 90, and so on
%
% For convenience, if time is a row vector and tilt/azimuth are column
% vectors (or vice versa), the results for GI, DNI, and DI will be expanded
% to a full matrix over all the values of each.

% This function is now the standard way to get clear sky irradiance.  The older slrGetCSKgi is deprecated and will be removed at some point in the indefinite future.

%% Find clear sky irradiance
% Check the args
if(nargin < 3 || isempty(tilt) )
	tilt = 0;
end
if(nargin < 4 || isempty(azimuth) )
	azimuth = 180;
end

% Save the time, then model GI
csk.time = time;
[csk.ghi, csk.sun] = getGlobalClearSkyIrradiance(csk.time, pos);
if(numel(tilt)==1 || isequal(size(tilt), size(time)))
	[csk.gi, csk.dni, csk.di] = ghi2gi( azimuth , tilt , pos , csk.time , csk.ghi );
elseif(all(min([size(tilt);size(time)])==[1 1])) % one row vector and one column vector
	st = size(time); s = size(tilt);
	[csk.gi, csk.dni, csk.di] = ghi2gi( repmat(azimuth,st) , repmat(tilt,st) , pos , repmat(csk.time,s) , repmat(csk.ghi,s) );
end

end
