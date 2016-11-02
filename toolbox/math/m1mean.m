function m = m1mean(m,dim)
% m1mean works like nanmean, except for integers:
%	for signed integers, -1 indicates nan
%	for unsigned integers, 0 indicates nan
%	double precision floating points are handled in the same way as signed integers

% generate the mask depending on the class of the input
if(regexp(class(m),'^uint\d+$'))
	% don't need to remove 'nan's in this case because they are by definition 0
	mask = m~=0;
else
	mask = m~=-1;
	% remove all the 'nan' values from the calculation
	m(~mask) = 0;
end

if(nargin == 2) %explicit dimension
	m = sum(m,dim);
	num = sum(mask,dim);
else % let sum() decide dimension
	m = sum(m);
	num = sum(mask);
end

% divide out:
m = m./num;

end
