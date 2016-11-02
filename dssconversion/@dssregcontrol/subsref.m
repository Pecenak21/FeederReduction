function [v varargout] = subsref(s,index)
% subsref allows us to look at object properties using . and ()
% subscripting

% Our strategy is just to follow the subscripts down the reference chain
% iteratively until we get to something that's not our class, at which
% point we hand off to that object.

% because I decided to use v instead of s for the value...
v = s;

for i = 1:length(index)
	% if we've gotten to the point where it's another class, just pass it
	% on to the next step
	if(~isa(v,mfilename('class')))
		v = subsref(v,index(i:end));
		break;
	end
	switch(index(i).type)
		case '()' %array indexing, just look up the nth item in the array
			v = v(index(i).subs{:});
		case '.'
			% struct indexing.  use 'get' for our class
			v = get(v,index(i).subs);
	end
end

% sometimes subscripting allows us to generate many outputs, which can be
% fed into a cell array or assigned to many variables, or whatever.  Here
% we handle that:
if(nargout>1)
	varargout = v(2:end);
	v = v{1};
end

end
