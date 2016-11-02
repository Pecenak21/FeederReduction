%% Matlab Library
%  Bryan's java library
%
%  Title: Time Linspace
%
%  Author: Bryan Urquhart
%
%  Description:
%    Takes two bu.util.Time objects and constructs a linear vector out of it,
%    similar to the linspace function
%
function timeval = jBuTimeLinspace( t1 , t2 , varargin )
%% Process Input Arguments

if( ~isa( t1 , 'bu.util.Time' ) )
  error( 'jBuTimeLinspace: Arg 1 must be a bu.util.Time.' );
end
if( ~isa( t2 , 'bu.util.Time' ) )
  error( 'jBuTimeLinspace: Arg 2 must be a bu.util.Time.' );
end

interval = 1; %[s]

if( ~isempty( varargin ) )
  args = argHandler( varargin );
  for idx =1:size( args )
    switch( args{idx,1} )
      case 'interval'
        interval = args{idx,2};
    end
  end
end

interval = interval/86400;


%% Generate vector

t1 = datenum( bu.util.Time.timeToDatevec( t1 )' );
t2 = datenum( bu.util.Time.timeToDatevec( t2 )' );

timeval = bu.util.Time.datevecToTime( datevec( t1:interval:t2 ) );
