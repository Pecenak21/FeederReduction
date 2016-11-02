%% Matlab Library
%
%  Title: What exactly is the matrix
%
%  Author: Morphius
%
%  Description:
%    Determine whether or not Neo is in the matrix
%
function ismatrix = isMatrix( A )

ismatrix = true;
if( ndims(A) ~= 2 )
  ismatrix = false;
else
  [m n] = size(A);
  if( sum([m n]>[0 0]) ~=2 ), ismatrix = false; end
end