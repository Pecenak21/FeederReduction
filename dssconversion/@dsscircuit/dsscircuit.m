function obj = dsscircuit(s,varargin)
% OpenDSS circuit object is based directly on the vsource object
cn = mfilename('class');
if(nargin == 0)
obj = dssvsource();
elseif(nargin == 1)
	obj = dssvsource(s);
else
	obj = dssvsource(s,varargin{:});
end
% make it an object
obj = class(struct(obj),cn);

end