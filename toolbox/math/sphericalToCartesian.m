%% Matlab math library
%
%  Title: Cartesian Coordinate Conversion
%
%  Author: Bryan Urquhart, and ?? 
%          -> someone else changed this function so it no longer
%             returns unit vectors! I am not sure I am the author
%             actually...
%
%  Description:
%    Returns the cartesian vector mapping for a given spherical coordinate
%    input.
%
%  Note: This function uses the y axis as a reference for azimuth when it
%  is the x axis that should be the reference. This is for surface
%  cartography applications obviously but the fcuntion is deceptively named
%  because it cannot be used in general for spherical to cartesian mapping.
%  
function [coord] = sphericalToCartesian( radius , azimuth , zenith )

% not clear why the cos/sin was changed to include a subtraction instead of
% using the other trig function.
coord.x = radius .* cos( pi/2 - azimuth ) .* sin( zenith );
coord.y = radius .* sin( pi/2 - azimuth ) .* sin( zenith );
coord.z = radius .* cos( zenith  );