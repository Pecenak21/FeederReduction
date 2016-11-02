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
%   The model only works up to 86 zenith. Anything beyond that is erroneous.
%
function [gi dni di sun] = ghi2gi( azimuth , tilt , position , time , ghi , varargin )
%% Process Input Arguments

% Verify that time and ghi are the same length
if( length(time) ~= length( ghi ) )
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

gi = [];
dni = [];
di = [];

%ZENITHLIMIT = 84.3;
ZENITHMAX_DIR   = 86.75;
ZENITHMAX_DIF   = 88.5;
ALPHAMAX_DIR   = 86.75;
ALPHAMAX_DIF   = 88.5;
sun = [];

% Process varargin
if( ~isempty( varargin ) )
  args = argHandler( varargin );
  for idx = 1:size(args,1)
    switch( args{idx,1} )
      case 'sun'
        sun = args{idx,2};
    end
  end
end



%% Compute the solar angles

% Use the NREL solar position algorithm
if( isempty( sun ) )
  sun = siSunPosition( time , position );
end

% Pull out values
So = bu.science.astronomy.Earth.SOLAR_CONSTANT.double() ./ sun.earthsundistance;

%% Compute the panel-solar angle
% Similar to sun-pixel angle

alpha = sky_imager_sun_pixel_angled( tilt , azimuth , sun.zenith , sun.azimuth );
sun.alpha = alpha;
%disp(sun.alpha)

%% Compute the diffuse function
%  This method is the Muneer method as documented by J. Page in the Practical
%  Handbook of photovoltaics. Here we are ignoring the contributions of
%  reflected radiation and are using the model to transpose diffuse horizontal
%  to diffuse (tilted). The diffuse function provides part of the formulation
%  for the diffuse radiation on a tilted plane, with the other requiring the
%  diffuse fraction as generated from the Boland et al. model below.

% --------------------------------------------
% This value can be tweaked if cloud fraction is known
b = 1.68;

% Compute the terms in the diffuse function
f_t1 = cos( B / 2 ).^2;
f_t2 = 2 * b / ( pi * ( 3 + 2 * b ) );
f_t3 = sin( B ) - B .* cos( B ) - pi * sin( B / 2 ).^2 ;

% Combine the terms of the diffuse function.
f = f_t1 + f_t2 * f_t3;
% --------------------------------------------

%% Compute the diffuse fraction
%  The diffuse fraction is generated from the Boland model which builds upon
%  years of work in generating a value for diffuse radiation simply from a value
%  of global irradiation.

kc = ghi ./ ( So .* cosd( sun.zenith ) );

dt = 1 ./ ( 1 + exp( -5 + 8.6 * kc ) );

dhi = ghi .* dt;

dhi( sun.zenith >= ZENITHMAX_DIF ) = ghi( sun.zenith >= ZENITHMAX_DIF );
% dhi( alpha      >= ALPHAMAX_DIF ) = ghi( alpha      >= ALPHAMAX_DIF );

%% Compute the beam irradiance and transmissivity

% Beam irradiance
dni = ghi .* ( 1 - dt ) ./ cosd( sun.zenith );
dni( sun.zenith >= ZENITHMAX_DIR ) = 0;
% dni( alpha      >= ALPHAMAX_DIR ) = 0;

% Beam transmissivity
kb = dni ./ So;

%% Compute the ground reflected radiation

reflectivity = 0.20;
groundslopefactor = ( 1 - cos( B ) ) / 2;
rgh = reflectivity * groundslopefactor .* ghi;

%% Compute the diffuse and global components

di = ( ( f .* ( 1 - kb ) + kb .* cosd( alpha ) ./ cosd( sun.zenith ) ) ) .* dhi;

% % Compute for solar elevation angles above 5.7
% di( sun.zenith < ZENITHLIMIT ) = ( ( f * ( 1 - kb( sun.zenith < ZENITHLIMIT ) ) ...
%   + kb( sun.zenith < ZENITHLIMIT ) .* cosd( alpha( sun.zenith < ZENITHLIMIT ) ) ...
%   ./ cosd( sun.zenith( sun.zenith < ZENITHLIMIT ) ) ) ) .* dhi( sun.zenith < ZENITHLIMIT );
% % Compute for solar elevation angles below 5.7
% di( sun.zenith >= ZENITHLIMIT ) = cos( B / 2 ) * ...
%   ( 1 + kb( sun.zenith >= ZENITHLIMIT ) * sin( B / 2 )^3 ) .* ...
%   ( 1 + kb( sun.zenith >= ZENITHLIMIT ) .* ...
%   cosd( alpha( sun.zenith >= ZENITHLIMIT ) ).^2 ...
%   .* sind( sun.zenith( sun.zenith >= ZENITHLIMIT ) ).^3 ) ...
%   .* dhi( sun.zenith >= ZENITHLIMIT );
% di = di';

% gi = dni .* cosd( alpha ) + di + rgh;
gi = dni .* cosd( alpha ) .* double(alpha<90) + di + rgh;

