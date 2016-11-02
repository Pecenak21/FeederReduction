function phase = getBusPhase(bus)
% return the bus phase with '.'
id = strfind(bus,'.');
if ~isempty(id)
    phase = bus(id(1):end);
end
end