%% Matlab Library
%
%  Title: Diff/Grad computation
%
%  Author: Bryan Urquhart
%
%  Description:
%    Computes diff, but over an interval. Resulting time series is the same
%    length. Starting from the first point each index corresponds to the
%    difference specified.
%

function y = lmGrad( x , varargin )
%% Process Input Arguments

% Set the default interval for computation
interval = 1;

if( ~isempty( varargin ) )
  args = argHandler( varargin );
  for idx = 1:size( args , 1 )
    switch( args{idx,1} )
      case 'interval'
        interval = args{idx,2};
    end
  end
end

% === Error checks ===
% Ensure the number of points does not exceed the length of the data set
if( interval > length(x) )
  warning( 'lmGrad: computation interval is greater than the length of the provided data set.' );
end

%% Set up processing kernel

mfilt = [-1 zeros(1,interval-1) 1];

y = conv( -mfilt , x ) / interval;
y = [y(2:end-interval);mean(y(2:end-interval))];

