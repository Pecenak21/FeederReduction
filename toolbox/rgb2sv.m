function [s,v] = rgb2sv(r)

if isa(r, 'uint8'),
	r = double(r) / 255;
elseif isa(r, 'uint16')
	r = double(r) / 65535;
end
	 
v = max(r,[],3);
if( isa(r,'double') )
    v(isnan(v)) = 0;
end
s = (v - min(r,[],3))./v;
s(~(v)) = 0; % removes NaNs

end
