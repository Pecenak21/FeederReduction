%% Solar Resource Assessment
%  Solar Library
%
%  Title: Relative opitcal airmass
%  
%  Author: Bryan Urquhart
%
%
%  Description:
%    Computes the relative optical airmass using the Kasten and Young (1989)
%    formula.
%
%
function am = airmass( zenith )

ZENITHMAX = 90;

zenith( zenith > ZENITHMAX ) = ZENITHMAX;

elevation = 90 - zenith;

a = 0.50572;
b = 6.07995;
c = 1.6364;

am = 1 ./  ( sind( elevation ) +  a * ( elevation + b ).^ -c );