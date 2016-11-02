%% Matlab Library
%  BU Java Libraries
%
%  Title: Separate a time period into days
%
%  Author: Bryan Urquhart
%
%  Description:
%    Takes a time range and converts it into an array of single days
%

function [beginarr endarr ndays]  = buTimeDaySplit( t1 , t2 , varargin )
%% Process Input Arguments

% Store inputs into standard time struct
time.begin = t1;
time.end   = t2;
time.zone  = 'UTC';

if( ~isempty( varargin ) )
  args = argHandler( varargin );
  for idx = 1:size(args,1)
    switch( args{idx,1} )
      case 'timezone'
        time.zone = args{idx,2};
    end
  end
end


%% Split into days

% Determine the number of distinct days
t_ = time.end.subtract( time.begin );
unixtime = t_.unixTime(); clear t_;
ndays = ceil( unixtime / ( 3600 * 24 ) );

if( ndays == 0 )
  beginarr = t1;
  endarr = t2;
  ndays = 1;
  return;
end

% Allocate arrays to hold beginning and ending periods
t_.begin = javaArray('bu.util.Time',ndays);
t_.end   = javaArray('bu.util.Time',ndays);

% Construct array of 'begin' and 'end' times for individual day processing
for idx = 1:ndays
  switch( idx )
    
    case 1
      
      t_.begin(idx) = time.begin;
      t_.end  (idx) = t_.begin(idx).addDays(1).yearDOY().subSeconds(1);
      
    case ndays
      
      t_.begin(idx) = t_.begin(idx-1).addDays(1).yearDOY();
      t_.end  (idx) = time.end.subSeconds(1);
      
    otherwise
      
      t_.begin(idx) = t_.begin(idx-1).addDays(1).yearDOY();
      t_.end  (idx) = t_.begin(idx).addDays(1).yearDOY().subSeconds(1);
      
	end
	% Get the middle of the day to find if day is DST (US only)
  tMiddle = bu.util.Time( t_.begin(idx) ); tMiddle.hour = 12;
  t_.begin(idx) = t_.begin(idx).toUTC(bu.util.Time.toDST(tMiddle,time.zone));
  t_.end  (idx) = t_.end  (idx).toUTC(bu.util.Time.toDST(tMiddle,time.zone));
end

% Store in time struct
beginarr = t_.begin;
endarr   = t_.end;

% Clean up workspace
clear t_ unixtime;
