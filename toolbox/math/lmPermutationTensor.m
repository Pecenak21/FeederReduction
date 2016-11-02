%% Math library
%
% Title: Construct permutatuib tensor
%
% Author: gnovice
%         If something doesn't work, blame Bryan Urquhart
%
% Description:
%   Constructs an n-dimensional permutation tensor
%
function pTens = lmPermutationTensor( n )
[matgrid{1:n}] = ndgrid(1:n);
pairsidx = nchoosek(1:n,2);
pTens = sign(prod(cat(n+1,matgrid{pairsidx(:,2)})-cat(n+1,matgrid{pairsidx(:,1)}),n+1));