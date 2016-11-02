%% Matlab Library
%
%  Title: Convert timezone of bu time array to local or to UTC
%
%  Description:
%    Convert entire time array to timezone of choice
%
function time = buTimeZone( time , varargin );
%% Process Input Arguments

args = argHandler( varargin );

toUTC = [];
zone = '';
for idx = 1:size(args,1)
  switch( args{idx,1} )
    case 'change'
      switch( args{idx,2} )
        case 'local'  
          toUTC = false;
        case 'utc'
          toUTC = true;
      end
    case 'zone'
      zone = upper( args{idx,2} );
  end
end

if( isempty( toUTC ) )
  error( 'You must supply a direction of time shift (local/utc)' );
end
if( isempty( zone ) )
  error( 'You must provide a zone argument' );
end

%% Change times

for idx = 1:time.length
  if( toUTC )
    time(idx) = time(idx).toUTC(zone);
  else
    time(idx) = time(idx).toLocal(zone);
  end
end
