function [r, x, r1, x1, r2, x2] = zconv(r, x, r1, x1, r2, x2)
%zconv converts between matrix impedance and sequence impedances
%
% the exact operation performed depends on the number of inputs/outputs:
% To convert from sequence impedance to matrix impedance, do one of the
% following:
%		z = zconv(z0,z1);
%		z = zconv(z0,z1,z2);
%		z = zconv(r0,x0,r1,x1,...);
%		[x r] = zconv(z0,z1,...);
%		[x r] = zconv(r0,x0,r1,x1,...);
% If negative sequence impedance is not specified, it is assumed to be
% equal to positive sequence impedance.
%
% To convert from matrix impedance to sequence impedances, specify z or
% both r and x as inputs, and get r0,x0,r1,x1,r2,x2 or z0,z1,z2 as outputs.
% Keep in mind that not all impedance matrices correspond exactly to a set
% of sequence impedances.
a = (-1)^(-2/3);
ma = [1 1 1; 1 a a^2; 1 a^2 a];
if(numel(r)==1) % have series want matrix
	if(nargin==2) % complex inputs, zero and pos only
		z = diag([r, x, x]);		
	elseif(nargin==3) % complex inputs
		z = diag([r, x, r1]);		
	elseif(nargin==4) % separated inputs, zero and pos only
		z = diag([r+i*x, r1+i*x1, r1+i*x1]);
	elseif(nargin==6) % separated inputs zero, pos, neg
		z = diag([r+i*x, r1+i*x1, r2+i*x2]);
	end
	z = ma*z*ma'/3;
	if(nargout<=1)
		r = z;
	elseif(nargout==2)
		r = real(z);
		x = imag(z);
	else
		error('too many outputs');
	end
else %matrix to series
	if(nargin==1) % complex input
		z = r;
	elseif(nargin==2) % split inputs
		z = r + i*x;
	else
		error('invalid input. matrix type assumed, but that doesn''t seem right either');
	end
	z = ma'*z*ma/3;
	r = z([1 5 9]);
	if(abs(sum(z([2 3 4 6 7 8])))>abs(sum(r)*.1))
		warning('zconv:sketchyresult','Trying to convert a matrix impedance to series form when it doesn''t seem to fit that form!');
	end
	if(any(nargout==[2 4]) && (abs(r(2)-r(3))/abs(r(2))>1e-3))
		warning('zconv:sketchyresult','Requested output assumes z1=z2 but this is not the case');
	end
	if(nargout <= 3) % complex sequence impedances
		r=z(1); x=z(5); r1=z(9);
	elseif(nargout ==4 || nargout == 6) % split sequence impedances
		r2 = real(r(3)); x2 = imag(r(3));
		r1 = real(r(2)); x1 = imag(r(2));
		x = imag(r(1)); r = real(r(1));
	else
		error('invalid number of outputs');
	end
end

end
