%% convert from psudo coordinates to lat,lon coordinates
% Generate deployment configuration with ground map information for Fallbrook feeder 520
% load feeder520 circuit data. This data is generated from convert2dss tool
c = load('data/f909.mat','c'); c = c.c;
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
% bus/node: 09099435  -> real coords: (x,y) or (lon,lat): [-116.93067500000001°, 33.26428888888889°]
[x, y] = ismember('09099435',bl.id);
i = 1;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-116.93067500000001, 33.26428888888889];

% bus/node: 0909 -> real coords: (x,y) or (lon,lat): [ -117.01785833333334°, 33.23080277777778°]
[x, y] = ismember('0909',bl.id);
i = 2;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [ -117.01785833333334, 33.23080277777778];

% bus/node: 09093928 -> real coords: (x,y) or (lon,lat): [-116.98524166666667°,33.220394444444445°]
[x, y] = ismember('09093928',bl.id);
i = 3;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-116.98524166666667,33.220394444444445];

% bus/node: 0909103  -> real coords: (x,y) or (lon,lat): [-116.95702222222222°, 33.268880555555555°]
[x, y] = ismember('0909103',bl.id);
i = 4;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-116.95702222222222, 33.268880555555555];

% bus/node: 09097319 -> real coords: (x,y) or (lon,lat): [-116.95795555555556°, 33.23663055555556°]
[x, y] = ismember('09097319',bl.id);
i = 5;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord =[-116.95795555555556, 33.23663055555556];

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
buslist2_909(length(bl.coord)) = bus(1);
for i = 1:length(buslist2_909)
	buslist2_909(i).id = bl.id(i);
	buslist2_909(i).psudo_coord = bl.coord(i,:);
	buslist2_909(i).coord(1) = polyval2(p2lon,bl.coord(i,1),bl.coord(i,2));
	buslist2_909(i).coord(2) = polyval2(p2lat,bl.coord(i,1),bl.coord(i,2));
end
buslist2 = buslist2_909;
temp = 'tmp/f480';
bl2 = dsswrite(c,'buslist2_909',1,[temp 'buslist2_909']);
save('data/buslist2_909.mat','buslist2');
save('data/buslist2_' c.circuit.Name '.mat','buslist2');
%% convert from lat,lon to psudo 