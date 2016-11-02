%% Matlab Library
%
%  Title: Time Vector Index
%
%  Author: Bryan Urquhart
%
%  Description:
%    Get the index of a given time point in a vector.
%
function index = jIndexOfTime( timearr , time , varargin )
%% Process Input Args

exact = true;

if(~isempty(varargin))
  args = argHandler( varargin );
  for idx = 1:size(args,1)
    switch( args{idx,1} )
      case 'exact'
        exact = args{idx,2};
    end
  end
end

%% Search

if( exact )
  index = java.util.Arrays.binarySearch( timearr , time );
  return;
end


index = -1;
for idx = 1:timearr.length
  
  % We first check to see that the time is not less than the first index
  % and if it is we quit and return that the time wasn't found'
  if( idx == 1 )
    if( time.compareTo( timearr( 1 ) ) < 0 ), break; end;
    index = 1;
  end
  
  % Check the current time at idx is greater than (later than) the cutoff time.
  % If it is we haven't looked far back enough yet so we need to continue
  % looping.
  if( timearr(idx).compareTo( time ) < 0 ), continue; end
  
  % When we don't get caught by the above continue, we are far back enough in
  % the data and we can store the index as our period of retrospective analysis.
  index = idx;
  break;
end
