function [c,indA] = grouping_unique(a)
% This is a modified copy of Matlab R2012a's implementation of 'unique'
%  I needed a copy of unique that would return indA (the second return) as
%  a cell array containing all the indices that match each unique value,
%  rather than just one
% 'R2012a' flag implementaion

% Determine if A is a row vector.
rowvec = isrow(a);

% Convert to column
a = a(:);
numelA = numel(a);

% Sort A and get the indices needed.
[sortA,indSortA] = sort(a);

% groupsSortA indicates the location of non-matching entries.
if isnumeric(sortA) && (numelA > 1)
	dSortA = diff(sortA);
	if (isnan(dSortA(1)) || isnan(dSortA(numelA-1)))
		groupsSortA = sortA(1:numelA-1) ~= sortA(2:numelA);
	else
		groupsSortA = dSortA ~= 0;
	end
elseif(iscellstr(sortA))
	groupsSortA = ~strcmp(sortA(1:numelA-1), sortA(2:numelA));
else
	groupsSortA = sortA(1:numelA-1) ~= sortA(2:numelA);
end

if (numelA ~= 0)
	groupsSortA = [true; groupsSortA];          % First element is always a member of unique list.
else
	groupsSortA = zeros(0,1);
end

% Extract unique elements.
c = sortA(groupsSortA);         % Create unique list by indexing into sorted list.

% Find indA.
if nargout <= 1
	warning('customUnique:missingThePoint','you should be using the bulit-in unique function if you only need one output');
end
x = find(groupsSortA);
% This is the magic sauce of our function that matlab's leaves out.
% Instead of just taking indA=indSortA(x), we grab the full range for each
% unique value
y = [x(2:end)-1; numel(groupsSortA)];
indA = cell(size(x));
for i=1:length(x)
	indA{i} = indSortA(x(i):y(i));
end


% If A is row vector, return C as row vector.
if rowvec
	c = c.';
end

end
