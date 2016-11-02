%% Matlab Library
%  Bryan's java library
%
%  Title: Position conversion
%
%  Author: Bryan Urquhart
%
%  Description:
%    Takes a struct with longitude, latitude, and altitude fields and
%    returns a bu.science.geography.Position object
%
function pOut = jToPosition( pIn )

if( isa( pIn , 'bu.science.geography.Position' ) )
	pOut = pIn; % Caution: this does not copy the object
	return;
end

pOut = bu.science.geography.Position( pIn.longitude , pIn.latitude , pIn.altitude );
