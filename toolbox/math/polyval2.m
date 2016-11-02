function z = polyval2(p,x,y)
% POLYVAL2 evaluates a 2 dimensional polynomial
% 
% Z = POLYVAL(P,X,Y) returns the value of a 2D polynomial P evaluated at (X,Y).  For a polynomial of order N, P
% is a vector of length (N+1)*(N+2)/2 containing the polynomial coefficients in
% ascending powers:
%
%   P = [p00 p10 p01 p20 p11 p02 p30 p21 p12 p03...]
%
% e.g. For a 2nd order fit, polyval2.m evaluates the matrix equation:
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
% NOTE: this is not the same polyval2 you will find on the file exchange.
% That version insists on evaluation on a grid, which is undesirable here

% Check input
if( any(size(x)~=size(y)) )
	error('polyval2:InvalidXY', 'X and Y must be same size.');
end
n=sqrt(2*numel(p)+0.25)-1.5;
if( ~isvector(p) || mod(n,1))
    error('polyval2:InvalidP',...
            'P must be a vector of length (N+1)*(N+2)/2, where N is order.');
end

s = size(x);
x = x(:);
y = y(:);
p = p(:);

% Construct "Vandermonde" matrix.
V=zeros(numel(x),numel(p));
V(:,1) = 1;
ordercolumn=1;
for order = 1:n
    for ordercolumn=ordercolumn+(1:order)
        V(:,ordercolumn) = x.*V(:,ordercolumn-order);
    end
    ordercolumn=ordercolumn+1;
    V(:,ordercolumn) = y.*V(:,ordercolumn-order-1);
end

z=V*p;
z=reshape(z,s);
