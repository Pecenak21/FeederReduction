function [X index] = exciserows(X,method)
%% Remove any rows that contain NaN values

if nargin < 2
	if( iscell(X) )
		method = 'cell';
	else
		method = 'nan';
	end
else
	% Convert to lower case for comparison
	method = lower(method);
end

switch( method )
  case 'nan'
    if( nargout == 2 )
      index = any(isnan(X),2);
    end
    X(any(isnan(X),2),:) = [];
  case 'zero'
    if( nargout == 2 )
      index = (mean(X,2)==0);
    end
    X(mean(X,2)==0,:) = [];
  case 'cell'
    index = any(cellfun('size',X,1)==0,2);
    X(index,:) = [];
  case 'false'
    if( nargout == 2 )
      index = any(~X,2);
    end
    X(any(~X,2),:) = [];
  case 'true'
    if( nargout == 2 )
      index = any(X,2);
    end
    X(any(X,2),:) = [];
end

end