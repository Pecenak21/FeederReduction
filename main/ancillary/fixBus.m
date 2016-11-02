function [bus, nphase] = fixBus(c,bus)
% check to see if the bus is consistent with the circuit config and fix it if not.

% connect to the correct phase of the bus by finding the components that are connected to that bus
% get element(s) connected to that bus
if iscell(bus), bus = bus{1}; end
el = findElemByBus(c,bus,0);
if isempty(el), error('This bus is isolated. Need to remove it in feederSetup.m before running this.');end
if isfield(el,'load') % if there is a load then connect to this load
    % if there is a load that is divided in to different phases but is represented by multiple loads then connect PV to all phases of this load
    % check by comparing kws and kvars of the loads + no load is 3 phases
    if length(el.load)>1 && length(unique([el.load.Kw])) == 1 ...
            && length(unique([el.load.Kvar])) == 1 ...
            && sum([el.load.Phases]==3) == 0
        for j = 1:length(el.load)
            phase = getBusPhase(el.load(j).bus1);
            if isempty(strfind(bus,phase))
                bus = [bus phase];
            end
        end
        % count number of dots to figure out number of phases
        nphase = length(strfind(bus,'.'));
        if nphase == 0, nphase = 3; end
    else % else just connect to the bus of first load
        bus = el.load(1).bus1;
        nphase = el.load(1).phases;
    end
elseif isfield(el,'line')
    % just use the bus info of the first line
    if strfind(lower(el.line(1).bus1),lower(bus))
        bus = el.line(1).bus1;
    else
        bus = el.line(1).bus2;
    end
    nphase = el.line(1).phases;
else
    error('Not support this comp besides from load and line yet. Write code here to handle it.');
end

end