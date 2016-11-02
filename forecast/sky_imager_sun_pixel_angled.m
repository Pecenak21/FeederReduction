function alpha = sky_imager_sun_pixel_angled( zenith , azimuth , sol_zenith , sol_azimuth )
% SKY_IMAGER_SUN_PIXEL_ANGLED computes the sun pixel angle mapping, in degrees
%
% Sun pixel angle is the angle between the sun and a given pixel.  It is
% just the angle between two vectors specified in spherical coordinates, so
% this function can be abused to calculate the angle between vectors that
% are not necessarily the sun and a pixel direction.
%
% Examples:
%       alpha = SKY_IMAGER_SUN_PIXEL_ANGLED( zenith , azimuth , sun.zenith, sun.azimuth )
%       alpha = SKY_IMAGER_SUN_PIXEL_ANGLED( zenith , azimuth , sun )
%
% See Also: siImager.sunPixelAngles

%% Check Inputs
if(nargin == 2)
	sol_zenith = azimuth.zenith;
	sol_azimuth = azimuth.azimuth;
	azimuth = zenith.azimuth;
	zenith = zenith.zenith;
elseif(nargin == 3)
	sol_azimuth = sol_zenith.azimuth;
	sol_zenith = sol_zenith.zenith;
end

%% Compute using spherical cosine formula
% Thanks to Shahrouz Alimo for suggesting this formula
cos_alpha = cosd(sol_zenith) .* cosd(zenith) + sind(sol_zenith) .* sind(zenith) .* cosd(azimuth - sol_azimuth);
% Take the arc cosine and get the result in degrees
alpha = acosd( cos_alpha );
