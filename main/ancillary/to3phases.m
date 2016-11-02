function c = to3phases(c)

for i = 1:length(c.load)
	c.load(i).Phases = 3;
	c.load(i).bus1 = cleanBus(c.load(i).bus1);
	c.load(i).kv = 12;
end

% %% line: to 3 phases
for i = 1:length(c.line)
	c.line(i).Phases = 3;
	c.line(i).bus1 = cleanBus(c.line(i).bus1);
	c.line(i).bus2 = cleanBus(c.line(i).bus2);
end

% generator
if isfield(c,'generator')
    for i = 1:length(c.generator)
        c.generator(i).Phases = 3;
        c.generator(i).bus1 = cleanBus(c.generator(i).bus1);
    end
end

% pvsystem
if isfield(c,'pvsystem')
    for i = 1:length(c.pvsystem)
        c.pvsystem(i).Phases = 3;
        c.pvsystem(i).bus1 = cleanBus(c.pvsystem(i).bus1);
    end
end

% transformer
if isfield(c,'transformer')
    for i = 1:length(c.transformer)
        c.transformer(i).Buses = cleanBus(c.transformer(i).Buses);
        c.transformer(i).Phases = 3;
    end
end

end