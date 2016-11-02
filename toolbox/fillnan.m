function y = fillnan(x)
% fill nan values in array with most recent value in the array

y = x;
% nan values id
a = find(isnan(x));

% not nan values id
na = find(~isnan(x));

% patch the nan id with previous not nan values
for i = 1:length(a)
	id = find( na < a(i), 1, 'last' );
	if ~isempty(id)
		y(a(i)) = x(na(id));
	end
end
