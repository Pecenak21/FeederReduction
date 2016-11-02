% Math library
%
% Title: Construct rotation tensor
%
% Author: Bryan Urquhart
%
% Description:
%   Using the permuation tensor and a matrix this code constructes a
%   rotation tensor from the input vector.
%
%% Function header
function R = lmRotationTensor( x )
%% Check for symbolic var

if( isa(x,'sym') )
	y = x;
	x = [1;2;3];
	
end

n = length(x);
pTens = lmPermutationTensor( n );
X = repmat( reshape( x , [ones(1,n-1) n]) , [n*ones(1,n-1) 1] );
R = -sum( pTens .* X , n);

if( exist('y','var' ) )
	s = sign(R); % Get the signs
	R = abs(R); % Convert indices to positive numbers
	z = R == 0;
	R = y( R + z ); % Index into the sym vector, handle zeros
	R = s.*R.*~z; % Change signs and zero out entries
end
