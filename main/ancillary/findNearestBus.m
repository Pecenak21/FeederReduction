function [bus, nphase] = findNearestBus(c,bus)
% this function works for 5 SDG&E feeders only. Output includes specifying phases.

bus = lower(cleanBus(bus));
%% Load buslist2 for the SDG&E c
c_name =  c.circuit.Name;
try
    buslist2 = load(['data/buslist2_' c_name '.mat']);
catch
    buslist2 = load([normalizePath('$KLEISSLLAB24-1')  '/database/gridIntegration/PVimpactPaper/buslist/buslist2_' c_name '.mat']);
end
buslist2 = buslist2.buslist2;
% find bus id
try 
    id = ismember(lower({buslist2.id}'),bus);
catch 
    id = ismember(lower([buslist2.id]'),bus);
end
loc_geo = buslist2(id).coord;
%% only use buses in buslist2 that exist in the c
try 
    id = ismember(lower({buslist2.id}'),lower(c.buslist.id));
catch 
    id = ismember(lower([buslist2.id]'),lower(c.buslist.id));
end
buslist2 = buslist2(id);
%% Find closest bus location based on geodistance.
d = reshape([buslist2.coord],[2 length(buslist2)])';
dist = (d(:,1)-loc_geo(1)).^2 + (d(:,2)-loc_geo(2)).^2;

[~, i] = min(dist);
bus = buslist2(i).id;
[bus, nphase] = fixBus(c,bus);
end