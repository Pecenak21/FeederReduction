function d = insertNaNs(d, i)
% INSERTNANS inserts NaN values at specified indices
% This is most likely useful to insert breaks into a line in a figure.  Multiple lines could be used, but that tends to be very slow compared to sticking in NaNs.
%
% When used with a matrix (as opposed to a vector), inserts rows of NaNs.
%
% Note that i identifies the indexes in the array as it currently is, rather than indexes after some nans are inserted.

% note whether we had a row vector or a column vector
wasrow = false;
if(isrow(d))
	wasrow = true;
	d = d(:);
end

s = size(d,2);

% convert to cells:
c = mat2cell(d,diff([1;i(:);size(d,1)+1]))';

% add nans
c(2,:) = {nan(1,s)};

% back to a matrix
d = vertcat(c{1:end-1}); % leave off the trailing nan

% if it was a row, switch it back!
if(wasrow)
	d = d';
end

end
