%% convert from psudo coordinates to lat,lon coordinates
% Generate deployment configuration with ground map information for Fallbrook feeder 520
% load feeder520 circuit data. This data is generated from convert2dss tool
c = load('data/f480.mat','c'); c = c.c;
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
% bus/node: 04809924  -> real coords: (x,y) or (lon,lat): [-117.2446833°, 32.7421833°]
[x, y] = ismember('04809924',bl.id);
i = 1;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-117.2446833, 32.7421833];

% bus/node: 04804780 -> real coords: (x,y) or (lon,lat): [ -117.2564472°, 32.7261500°]
[x, y] = ismember('04804780',bl.id);
i = 2;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-117.2564472, 32.7261500];

% bus/node: 04801426 -> real coords: (x,y) or (lon,lat): [-117.2495222°,32.7082167°]
[x, y] = ismember('04801426',bl.id);
i = 3;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-117.2495222,32.7082167];

% bus/node: 04803609  -> real coords: (x,y) or (lon,lat): [-117.2409222°, 32.7250861°]
[x, y] = ismember('04803609',bl.id);
i = 4;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-117.2409222, 32.7250861];

% bus/node: 04808232 -> real coords: (x,y) or (lon,lat): [-117.2493417°, 32.7303222°]
[x, y] = ismember('04808232',bl.id);
i = 5;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord =[-117.2493417, 32.7303222];

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
buslist2(length(bl.coord)) = bus(1);
for i = 1:length(buslist2)
	buslist2(i).id = bl.id(i);
	buslist2(i).psudo_coord = bl.coord(i,:);
	buslist2(i).coord(1) = polyval2(p2lon,bl.coord(i,1),bl.coord(i,2));
	buslist2(i).coord(2) = polyval2(p2lat,bl.coord(i,1),bl.coord(i,2));
end
temp = 'tmp/f480';
bl2 = dsswrite(c,'buslist2',1,[temp 'buslist2']);
save('data/buslist2_480.mat','buslist2');
%% convert from lat,lon to psudo 