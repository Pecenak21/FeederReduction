%% Matlab Library
%
%  Title: Diff/Grad computation
%
%  Author: Bryan Urquhart
%
%  Description:
%    Computes diff for an interval of length n
%

function y = lmDiff( x , n , varargin )
%% Process Input Arguments

% Set the default interval for computation
if( nargin == 1 )
	n = 1;
end

if( ~isempty( varargin ) )
  args = argHandler( varargin );
  for idx = 1:size( args , 1 )
    switch( args{idx,1} )
    end
  end
end

% === Error checks ===
% Ensure the number of points does not exceed the length of the data set
if( n > length(x) )
  warning( 'lmGrad: computation interval is greater than the length of the provided data set.' );
end

%% Set up processing kernel

mfilt = [-1 zeros(1,n-1) 1];

y = conv( -mfilt , x );
y = [y(2:end-n);mean(y(2:end-n))];

