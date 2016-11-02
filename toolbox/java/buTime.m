%% Matlab Library
%  BU Java Libraries
%
%  Title: Convert time into bu time
%
%  Author: Bryan Urquhart
%
%  Description:
%    Converts a matlab time format into a bu time format
%

function timeout = buTime( timein )
%% Process Input Arguments

% Process variations of bu time input arguments
if( isa( timein , 'bu.util.Time' ) )
  timeout = timein;
  return;
end

if( isa( timein , 'bu.util.Time[]' ) )
  timeout = timein;
  return;
end

classStr = java.lang.String( class( timein ) );
if( classStr.indexOf( 'bu.util.Time' ) >= 0 )
  timeout = timein;
  return;
end

%% Check for Datenum

if( isa( timein , 'double' ) && size( timein , 2 ) == 1 )
  dvec = datevec( timein );
  timeout = bu.util.Time.datevecToTime( dvec );
  return;
end

%% Check for datevec

if( ( isa( timein , 'double' ) || isinteger( timein ) ) && size( timein , 2 ) == 6 )
  timeout = bu.util.Time.datevecToTime( timein );
  return;
end

%% Unknown time type

if( isempty( timein ) )
  timeout = [];
  return;
end

error( 'Unknown time format specified' );

