%% Matlab Libary
%
%  Title: Root Mean Square Error
%
%  Author: Bryan Urquhart
%
%  Description:
%    Computes the element by element rms error with respect to zero or the
%    second parameter if one is given.
%
%
function [err] = rmse( x , y )
%% Process Input Arguments

if( nargin == 2 )
	if( size(x) ~= size(y) )
		error( 'Both arguments must have identical size' );
	end
	% Compute RMSE
	err = sqrt( nanmean( ( x(:) - y(:) ).^2 ) );
else
	% Compute RMSE
	err = sqrt( nanmean( ( x(:) ).^2 ) );
end

end
