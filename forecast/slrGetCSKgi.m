%% Solar Resource Assessment Library
% 
% Title: Clear sky GI structure
%
% Author: Bryan Urquhart
%
% Description:
%   Get complete clear sky structure including plane of array irradiance, dni,
%   diffuse, and all the solar angles.
%
%	=====Deprecated=====
%	This function has been deprecated in favor of clearSkyIrradiance() and will be removed at a future date
function csk =  slrGetCSKgi( pos , time , timezone , tilt , azimuth )
%% Proces Input Arguments

if(nargin < 3 || isempty(timezone) )
	timezone = 'UTC';
end
if(nargin < 4 || isempty(tilt) )
	tilt = 0;
end
if(nargin < 5 || isempty(azimuth) )
	azimuth = 180;
end


if( isa(time,'double') )
	time = bu.util.Time.datevecToTime(datevec(time(:)));
end

csk.timezone = timezone;
if( exist('timezone','var') && ~strcmp(timezone,'UTC') )
  for i = 1:numel(time)
    time(i) = time(i).toUTC(timezone);
  end
end

if( isempty( time ) )
  csk = [];
  return;
end

%% === GLOBAL HORIZONTAL IRRADIANCE ===

% User message
% fprintf( '\t1. Computing clear sky global horizontal irradiance ... ');
% timer_ = tic;

% TEMP---
% pos = position;
% time = timevec;

% Perform Computation
geo = geoGet();
if(isstruct(pos))
	pos = bu.science.geography.Position( pos.longitude, pos.latitude, pos.altitude );
end
if(numel(time)==1)
	csk.ghi = geo.csk.getGlobalClearSkyIrradiance( pos.longitude , pos.latitude , pos.altitude , time );
	sunfield = {};
	csk.time = time;
else
	csktmp = geo.csk.getGlobalClearSkyIrradiance( pos , time );
	% Copy from java to matlab format
	csk.time									= csktmp.time;
	csk.ghi										= csktmp.ghi;
	csk.sun.zenith						= csktmp.sun.zenith;
	csk.sun.azimuth						= csktmp.sun.azimuth;
	csk.sun.earthsundistance	= csktmp.sun.earthSunRadius;
	csk.sun.zenith						= csktmp.sun.zenith;
	sunfield = {'sun',csk.sun};
end


% Store solar fields for reuse
% These seem to be broken? but ghi2gi thinks it can calc them on its own, so maybe we'll let it? 
% Response from Bryan: it was broken. I modified the way this function works to be a bit faster (hopefully)
%csk.sun.time              = dat.csk.sSun.time;
%csk.sun.zenith            = dat.csk.sSun.zenith;
%csk.sun.azimuth           = dat.csk.sSun.azimuth;
%csk.sun.earthsundistance  = dat.csk.sSun.earthsundistance;

% User message
% fprintf( '%.2f s\n' , toc( timer_ ) );

%% === TRANSPOSITION TO GLOBAL IRRADIANCE ===
% User message
% fprintf( '\t2. Transposing to clear sky global irradiance ......... ' );
% timer_ = tic;

% Perform Computation
[csk.gi,csk.dni,csk.di,csk.sun] = ghi2gi( azimuth , tilt , pos , time , csk.ghi , sunfield{:} );

% User message
% fprintf( '%.2f s\n' , toc( timer_ ) );
