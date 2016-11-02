%% Sky Imager Library
%
%  Title: Cartesian Coordinate Conversion
%
%  Author: Bryan Urquhart
%
%  Description:
%    Returns the cartesian unit vector mapping for a given spherical coordinate
%    input
%  
function [coord] = sphericalToCartesiand( radius , azimuth , zenith )

coord.x = radius .* cosd( 90 - azimuth ) .* sind( zenith );
coord.y = radius .* sind( 90 - azimuth ) .* sind( zenith );
coord.z = radius .* cosd( zenith  );