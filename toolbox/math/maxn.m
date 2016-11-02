%% Matlab Math Library
% Title : Find the max in a N-D matrix
% Author: Bryan Urquhart
% Description:
%   Finds the indicial location of the max value in the input matrix
%
function [C maxi] = maxn( A )
%% Process Input Arguments
if( isvector(A) ),
  [C maxi] = max(A);
  return;
end

n = ndims(A);

%% Get the max

maxi = nan(n,1);

[B,i] = max(A,[],n);
[~,index] = maxn(B);

% Compute the indicial location of the max
c = num2cell( index );
ind = sub2ind(size(i),c{:});

maxi(1:end-1) = index;
maxi(end) = i(ind);

c = num2cell( maxi );
ind = sub2ind(size(A),c{:});  
C = A(ind);
