%% Solar Resource Assessment
%  Solar Radiation Library
%
%  Title: Horizontal to Tilted Transposition
%  
%  Author: Bryan Urquhart
%
%  Description:
%    This function will convert a time series of GHI and time to the Global,
%    direct and diffuse components on a horizontal surface. It uses the Muneer
%   (aka Page) model to perform the transposition, and it uses the Bowland
%   equation to generate the diffuse fraction from the global irradiance. The
%   conversion from tilted to horizontal is performed by the
%   complementary function
%
function [ghi dni dhi] = gi2ghi( azimuth , tilt , position , time , gi )
%% Process Input Arguments

% Verify that time and ghi are the same length
if( length(time) ~= length( gi ) )
  error( 'The input time and ghi must be the same length.' );
end

% Convert tilt to radians
B = tilt * pi / 180;

% Fix up the time
if(~isnumeric(time))
	if(ischar(time))
		time = datenum(time);
	else % assume some form of buTime
		% unfortunately, timeToDatevec returns the correct form for vector time and the transpose for scalar time, so:
		time = bu.util.Time.timeToDatevec(time);
		if(size(time,2)~=6), time = time'; end
		time = datenum(time);
	end
end

ghi = [];
dni = [];
dhi = [];

ZENITHLIMIT = 84.3;
ZENITHMAX   = 86;

%% Compute the solar angles

% Use the NREL solar position algorithm
sun = siSunPosition( time , position );

% Pull out values
So = bu.science.astronomy.Earth.SOLAR_CONSTANT.double() ./ sun.earthsundistance;
% for idx = 1:length(sunangles)
%   So(idx,1)      = bu.science.astronomy.Earth.SOLAR_CONSTANT.double() / sunangles(idx).earthsundistance;
%   sun.zenith(idx,1)  = sunangles(idx).zenith;
%   sun.azimuth(idx,1) = sunangles(idx).azimuth;
% end

%% Compute the panel solar zenith

alpha = sky_imager_sun_pixel_angled( tilt , azimuth , sun.zenith , sun.azimuth );

%% Compute the diffuse function
%  This method is the Muneer method as documented by J. Page in the Practical
%  Handbook of photovoltaics. Here we are ignoring the contributions of
%  reflected radiation and are using the model to transpose diffuse horizontal
%  to diffuse (tilted). The diffuse function provides part of the formulation
%  for the diffuse radiation on a tilted plane, with the other requiring the
%  diffuse fraction as generated from the Boland et al. model below.

% This value can be tweaked if cloud fraction is known
b = 1.68;

% Compute the terms in the diffuse function
f_t1 = cos( B / 2 )^2;
f_t2 = 2 * b / ( pi * ( 3 + 2 * b ) );
f_t3 = sin( B ) - B * cos( B ) - pi * sin( B / 2 )^2 ;

% Combine the terms of the diffuse function.
f = f_t1 + f_t2 * f_t3;

%% Compute the diffuse fraction
%  The diffuse fraction is generated from the Boland model which builds upon
%  years of work in generating a value for diffuse radiation simply from a value
%  of global irradiation.

kc = gi ./ ( So .* cosd( sun.zenith ) );

dt = 1 ./ ( 1 + exp( -5 + 8.6 * kc ) );

di = gi .* dt;

di( alpha      >= ZENITHMAX ) = gi( alpha      >= ZENITHMAX );
di( sun.zenith >= ZENITHMAX ) = gi( sun.zenith >= ZENITHMAX );

%% Compute the beam irradiance and transmissivity

% Beam irradiance
dni = gi .* ( 1 - dt ) ./ cosd( alpha );
dni( sun.zenith >= ZENITHMAX ) = 0;
dni( alpha      >= ZENITHMAX ) = 0;

% Beam transmissivity
kb = dni ./ So;

%% Compute the diffuse and global components

dhi = ( 1 ./ ( f * ( 1 - kb )  + kb .* cosd( alpha ) ./ cosd( sun.zenith ) ) ) .* di;

ghi = dni .* cosd( sun.zenith ) + dhi;
