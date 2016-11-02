function [ghi, sun] = getGlobalClearSkyIrradiance(time, latitude, longitude, altitude)
%getGlobalClearSkyIrradiance(time, latitude, longitude, altitude)
% computes

if(isstruct(latitude)) % user gave a position struct instead of lat/lon/alt
	pos = latitude;
	altitude = pos.altitude;
elseif(isa(latitude,'bu.science.geography.Position'))
    pos = struct(latitude);
    altitude = pos.altitude;
else
	pos = struct('latitude', latitude, 'longitude', longitude, 'altitude', altitude);
end
sun = siSunPosition(time, pos);
zenith = reshape(vertcat(sun.zenith),size(time));
esd = reshape(vertcat(sun.earthsundistance),size(time));

% Get the airmass
am = airmass( zenith );

% Compute the extraterrestrial solar constant for this time
solarConstant = 1366.1 ./ ( esd .^ 2 ); %bu.science.astronomy.Earth.SOLAR_CONSTANT has the value 1366.1

% Compute some alitude parameters included in the determination of
% clear sky transmissivity
fh1 = exp( - altitude / 8000.0 );
fh2 = exp( - altitude / 1250.0 );
cg1 = 0.0000509 * altitude + 0.868; %0.886;
cg2 = 0.0000392 * altitude + 0.0387;

% Aerosol Accounting
TL  = LinkeTurbidity( time, pos );

transmissivity = cg1 * exp( - cg2 * am .* (fh1+fh2*(TL-1)) ) .* exp( 0.01 * am .^ 1.8 );

ghi = transmissivity .* solarConstant .* cosd( zenith );
ghi(ghi<0) = 0;

end
