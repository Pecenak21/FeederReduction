%% Matlab Libary
%
%  Title: Mean Bias Error
%
%  Author: Bryan Urquhart
%
%  Description:
%    Computes the element by element mean bias error with respect to zero or the
%    second parameter if one is given.
%
%
function [err] = mbe( x , varargin )
%% Process Input Arguments

if( ~isempty( varargin ) )
  y = varargin{1};
  
  if( size(x) ~= size(y) )
    error( 'Both arguments must have identical size' );
  end
else
  y = zeros( size(x) );
end

%% Compute RMSE

err = nanmean( ( x(:) - y(:) ) );