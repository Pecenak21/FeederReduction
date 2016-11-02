function [ c ] = gen2pv( c )
% replace generator by PV system with power output curve as input
% and PV rated output = generator active power rating

if isa(c, 'dssgenerator')
    gen = c;
elseif isstruct(c)
    try gen = c.generator;
    catch e
        error('No pvsystem in input circuit!');
    end
else
    error('Wrong data type in!');
end

pv(length(gen)) = dsspvsystem();
for i = 1:length(gen)
	pv(i).kVA = gen(i).kw*1.1; 
    pv(i).pmpp = gen(i).kw;
    pv(i).irradiance = 1;
    pv(i).pf = 1;
    pv(i).Name = gen(i).Name;
	pv(i).kv = gen(i).kv;
	pv(i).bus1 = gen(i).bus1;
	pv(i).phases = gen(i).phases;
	% NO VAR SUPPORT
%  	gen(i).Kvar = 0;%.2*pv(i).kVA;
% 	gen(i).model = 7;
    pv(i).Enabled = gen(i).Enabled;
    pv(i).Daily = gen(i).Daily;
    pv(i).Yearly = gen(i).Yearly;
end

if isa(c,'dsspvsystem')
    c = pv;
else
    c.pvsystem = pv;
    c = rmfield(c,'generator');
end

end

