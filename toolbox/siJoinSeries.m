function s = siJoinSeries(newSeries,s)
% siJoinSeries joins different series files together 
%
% Inputs:
%	newSeries - the new element to be appended
%	s     - the series to appdend to (optional).  If this input is empty or unspecified, a new series is created.

% if s doesn't have anything in it, we return the original element, but with some things rotated or wrapped in cell arrays
if(nargin < 2 || isempty(s) )
	s = newSeries;
	return;
end

% otherwise we append the subitems one at a time
for fn = fieldnames(newSeries)'; fn = fn{1};
	s.(fn) = [s.(fn); newSeries.(fn)];
end

end