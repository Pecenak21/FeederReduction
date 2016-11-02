%% Math Library
%
%  Title: transpose 2D sub matrices of 3D array
%
%  Author: Bryan Urquhart
%
%  Description:
%    Performs transpose on 2D submatrices of a 3D array. This may be able
%    to be generalized, but I don't want to think about it now because I
%    only need this functionality. If someone sees this method and wants to
%    construct an N-D transpose with some sort of behavior I recommend the
%    name transposeN or transposen.

% Comment:
% I think that this all ends up being equivalent to
%
% B = permute(A, [2 1 3]);
%
% and that you could write the N-D version as
%
% B = permute(A, [2 1 3:ndims(A)]);
%
% Cheers, --Ben

%% Function Header
function B = transpose3( A )
%% Process input arguments

% Check number of dimensions
if( ndims(A) ~= 3 )
	error( 'Input matrix must have ndims(A) == 3.' );
end

% Get size
d = size(A);

%% Perform transpose

% Swap indices for transposition
tmp = d(1);
d(1) = d(2);
d(2) = tmp; clear tmp;

% Allocate output
B = nan(d);

% Transpose
for i = 1:d(3)
	B(:,:,i) = A(:,:,i)';
end
