function [busId] = locateBus(circuit, x_coordinate, y_coordinate)
% locate the bus that a component should be connected to based on its
% location.
%
% input:
%           circuit: opendss matlab struct circuit
%           x_coordinate
%           y_coordinate
% output:
%           bus : without specifying phases

if ischar(x_coordinate)
    x_coordinate = str2double(x_coordinate);
end

if ischar(y_coordinate)
    y_coordinate = str2double(y_coordinate);
end

%% Load buslist2 for the good circuit
circuit_name =  circuit.circuit.Name;
try
    buslist2 = load(['data/buslist2_' circuit_name '.mat']);
catch
    buslist2 = load([normalizePath('$KLEISSLLAB24-1')  '/database/gridIntegration/PVimpactPaper/buslist/buslist2_' circuit_name '.mat']);
end
buslist2 = buslist2.buslist2;
%% only use buses in buslist2 that exist in the circuit
try 
    id = ismember(lower({buslist2.id}'),lower(circuit.buslist.id));
catch
    id = ismember(lower([buslist2.id]'),lower(circuit.buslist.id));
end
buslist2 = buslist2(id);
%% Find closest bus location based on geodistance.
d = reshape([buslist2.coord],[2 length(buslist2)])';
d(:,1) = d(:,1) - x_coordinate;
d(:,2) = d(:,2) - y_coordinate;
d(:,3) = d(:,1).^2 + d(:,2).^2;

[~, i] = min(d(:,3));
busId = buslist2(i).id;

end