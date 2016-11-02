function p = polyfit2(x,y,z,n,w)
% POLYVAL2 fits a surface to a 2 dimensional polynomial
% 
% P = POLYFIT2(X,Y,Z,N,W) finds the (N+1)*(N+2)/2 of an Nth order 2
% dimensional polynomial that fits the data best in a least-squares sense
%
% P is a row vector starting with the lowest order coefficients, e.g.
%
%   P = [p00 p10 p01 p20 p11 p02 p30 p21 p12 p03...]
%
% e.g. For a 2nd order fit, polyfit2.m evaluates the matrix equation:
%
%    Z = V*P    or
%
%    Z = [1  x  y  x^2  xy  y^2]  [p00
%                                  p10
%                                  p01
%                                  p20
%                                  p11
%                                  p02]

%
%
% NOTE: this is not the same  as the polyfitweighted2 you will find on the file exchange.
% That version insists on evaluation on a grid, which is undesirable here

% e.g. For a 3rd order fit, 
% the regression problem is formulated in matrix format as:
%
%   wZ = V*P    or
%
%                      2       2   3   2     2      3
%   wZ = [w  wx  wy  wx  xy  wy  wx  wx y  wx y   wy ]  [p00
%                                                        p10
%                                                        p01
%                                                        p20
%                                                        p11
%                                                        p02
%                                                        p30
%                                                        p21
%                                                        p12
%                                                        p03]
%

x = x(:);
y = y(:);

if( any(size(x)~=size(y)) || any(size(y)~=size(z)) || (nargin > 4 && any(size(z)~=size(w))) )
    error('polyfit2:SizeMismatch','x,y,z,w must be of the same size')
end

x = x(:);
y = y(:);
z = z(:);

% Construct weighted Vandermonde matrix.
V=zeros(numel(z),(n+1)*(n+2)/2);
if(nargin > 4)
	w = w(:);
	V(:,1) = w;
	z = z.*w;
else
	V(:,1) = 1;
end
ordercolumn=1;
for order = 1:n
    for ordercolumn=ordercolumn+(1:order)
        V(:,ordercolumn) = x.*V(:,ordercolumn-order);
    end
    ordercolumn=ordercolumn+1;
    V(:,ordercolumn) = y.*V(:,ordercolumn-order-1);
end

% Solve least squares problem.
[Q,R] = qr(V,0);
ws = warning('off','all'); 
p = R\(Q'*z);    % Same as p = V\z;
warning(ws);
if size(R,2) > size(R,1)
	warning('polyfit2:PolyNotUnique', 'Polynomial is not unique; degree >= number of data points.')
elseif condest(R) > 1.0e10
	warning('polyfit2:RepeatedPointsOrRescale', 'Polynomial is badly conditioned. Remove repeated data points\n         or try centering and scaling as described in HELP POLYFIT.')
end
p = p.';          % Polynomial coefficients are row vectors by convention.
