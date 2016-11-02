function y = medfilt1mod(x, n, type, withweight)
% modified median filter 1D to do median filter 
%
% type :	'h' for history (not knowing the future, using history to get median value)
%			'c' centered median filter or 'normal' med filter

if ~exist('type','var') || isempty(type)
	type = 'c';
end

if ~exist('withweight','var') || isempty(withweight)
	ww = 0;
else
	ww = withweight;
end

nx = length(x);
if rem(n,2)~=1    % n even
    m = n/2;
else
    m = (n-1)/2;
end
if strcmp(type(1),'h')
	X = [nan(2*m,1); x(:)];
else
	X = [nan(m,1);x(:);nan(m,1)];
end
y = zeros(nx,1);

for i=1:nx
	xid = i + 2*m;
	xx = [];
	if ~ww
		xx = X(xid-n+1:xid);
	else
		for j = 1 : n 
			xx = [xx; ones(j,1)*X(xid -n + j)];
		end
	end
    y(i) = median(xx(~isnan(xx)));
end

%% replace the initial zero values by original values in the array (not doing med)
id = y(1:m*2)==0;
y(id) = x(id);
