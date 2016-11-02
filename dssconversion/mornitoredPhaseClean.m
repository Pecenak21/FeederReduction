function [p ps pfullform] = mornitoredPhaseClean(s)
% clean up mornitored phase.
% output:
%			p : {1,2,3} phase number; for delta/ll connection use the first of the two phases (so 1 for 1-2, 2 for 2-3, 3 for 3-1). Default to 1 if input is invalid.
%			ps: phase(s) monitored in form of 'xyz'
%			pfullform: mornitored phase represented as full form for OpenDSS bus. E.g: 1.2.3, 1.2, 2.3

if ischar(s)
	s(s==' ') = '';
	s = sort(lower(s));
	allphase = 'abcxyz';
	[val id] = ismember(s,allphase);
	ps = allphase(nonzeros(id));
	ps(ps=='a') = 'x';
	ps(ps=='b') = 'y';
	ps(ps=='c') = 'z';
	ps = unique(ps);
	switch ps
		case {'x','xy'}
			p = 1;
		case {'y','yz'}
			p = 2;
		case {'z','xz'}
			p = 3;
		case {'xyz'} % Using all three phases (p should be ignored)
			p = 0;
		otherwise
			warning('DSSConversion:mornitoredPhaseClean','Invalid monitored phase input');
			p = 1; ps ='x';
	end
elseif isnumeric(s)
	p = uint8(s);
	switch p
		case 1
			ps = 'x';
		case 2
			ps = 'y';
		case 3 
			ps = 'z';
		otherwise
			warning('DSSConversion:mornitoredPhaseClean','Invalid monitored phase input');
			p = 1; ps = 'x';
	end
else
	warning('DSSConversion:mornitoredPhaseClean','Invalid monitored phase input');
	p = 1; ps = 'x';
end
	pfullform = generateFullForm(ps);
end

function o = generateFullForm(ps)
	if length(ps) == 1, o = num2str(toNum(ps));
	elseif length(ps) > 1
		o = num2str(toNum(ps(1)));
		for i=2:length(ps)
			o = [o '.' num2str(toNum(ps(i)))];
		end
	else
		error('Invalid input');
	end
		
end

function o = toNum(p)
	switch p
		case {'a','x'}
			o = 1;
		case {'b','y'}
			o = 2;
		case {'c','z'}
			o = 3;
	end
end