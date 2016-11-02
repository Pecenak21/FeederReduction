function c = pv2gen(c)
% replace pvsystem by generator

if isa(c, 'dsspvsystem')
    pv = c;
elseif isstruct(c)
    try pv = c.pvsystem;
    catch e
        error('No pvsystem in input circuit!');
    end
else
    error('Wrong data type in!');
end

gen(length(pv)) = dssgenerator();
for i = 1:length(pv)
	gen(i).Kw = pv(i).kVA; 
	gen(i).Name = pv(i).Name;
	gen(i).kv = pv(i).kv;
	gen(i).bus1 = pv(i).bus1;
	gen(i).phases = pv(i).phases;
	% NO VAR SUPPORT
 	gen(i).Kvar = 0;%.2*pv(i).kVA;
	gen(i).model = 7;
    gen(i).Enabled = pv(i).Enabled;
    gen(i).Daily = pv(i).Daily;
    gen(i).Yearly = pv(i).Yearly;
end

if isa(c,'dsspvsystem')
    c = gen;
else
    c.generator = gen;
    c = rmfield(c,'pvsystem');
end
end