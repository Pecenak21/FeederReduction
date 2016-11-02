%% convert from psudo coordinates to lat,lon coordinates
% Generate deployment configuration with ground map information for Fallbrook feeder 520
% load feeder520 circuit data. This data is generated from convert2dss tool
c = load('data/f355.mat','c'); c = c.c;
%% we're only interested in bus location
bl = c.buslist;

% load from disk
% load('dssconversion/data/SDGE_feeders.mat');
% figure, hold on
% f = f520;
% for i = 1:length(f)
% 	plot([f(i).X]',[f(i).Y]');
% end

% let's manually match some buses from psudo-coords to real geographic coords
% do this using circuitvisualizer from convert2dss package and plotting
% function in SDGE_feeders.m
% bus/node: 03551325  -> real coords: (x,y) or (lon,lat): [-116.79521944444444°, 32.841613888888894°]
[x, y] = ismember('03551325',bl.id);
i = 1;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-116.79521944444444, 32.841613888888894];

% bus/node: 035560 -> real coords: (x,y) or (lon,lat): [ -116.77594444444445°, 32.81633055555555°]
[x, y] = ismember('035560',bl.id);
i = 2;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [ -116.77594444444445, 32.81633055555555];

% bus/node: 03552614 -> real coords: (x,y) or (lon,lat): [-116.76849722222222°,32.837808333333335°]
[x, y] = ismember('03552614',bl.id);
i = 3;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-116.76849722222222,32.837808333333335];

% bus/node: 0355506  -> real coords: (x,y) or (lon,lat): [-116.78528888888889°, 32.82876944444445°]
[x, y] = ismember('0355506',bl.id);
i = 4;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-116.78528888888889, 32.82876944444445];

% bus/node: 03554009 -> real coords: (x,y) or (lon,lat): [-116.7723°, 32.82743055555556°]
[x, y] = ismember('03554009',bl.id);
i = 5;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord =[-116.7723, 32.82743055555556];

%% 2D fitting
clear a b
for i = 1:length(bus)
	a(i).lon = bus(i).psudo_coord(1);
	a(i).lat = bus(i).psudo_coord(2);
	b(i).lon = bus(i).coord(1);
	b(i).lat = bus(i).coord(2);
end
%
p2lon = polyfit2([a.lon]',[a.lat]',[b.lon]',1);
p2lat = polyfit2([a.lon]',[a.lat]',[b.lat]',1);

% Generate other buses from the list
buslist2_355(length(bl.coord)) = bus(1);
for i = 1:length(buslist2_355)
	buslist2_355(i).id = bl.id(i);
	buslist2_355(i).psudo_coord = bl.coord(i,:);
	buslist2_355(i).coord(1) = polyval2(p2lon,bl.coord(i,1),bl.coord(i,2));
	buslist2_355(i).coord(2) = polyval2(p2lat,bl.coord(i,1),bl.coord(i,2));
end
buslist2 = buslist2_355;
temp = 'tmp/f355';
bl2 = dsswrite(c,'buslist2_355',1,[temp 'buslist2_355']);
save('data/buslist2_355.mat','buslist2');
save('data/buslist2_' c.circuit.Name '.mat','buslist2');
%% convert from lat,lon to psudo 