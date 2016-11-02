%% get circuit and buslist
c = feederSetup('pointloma','validated');
% we're only interested in bus location
bl = c.buslist;

%% get bus location in both psudo and lat-lon coordinates 
% psudo buslist location from SDGE data
load('SDGE_feeders.mat');
figure, hold on
f = f480;
for i = 1:length(f)
	plot([f(i).X]',[f(i).Y]');
end

% manually map psudo and real coordinates of a few buses 
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

% 2D fitting
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

%% Load PV data. The data load won't be in object types like in convert2dss package but we can deal with that.
pv2 = load([normalizePath('$KLEISSLLAB24-1')  '/database/gridIntegration/PVimpactPaper/virtualPV/f480scenario2_pv.mat']);
pv2 = pv2.pv;
pv = struct();
for i = 1:length(pv2)
    pv(i).bus1 = pv2(i).bus1;
    pv(i).kVA = pv2(i).kVA;
    pv(i).Name = pv2(i).Name;
end

% add useful fields
pv().lon_psudo = [];
pv().lon = [];
pv().lat_psudo = [];
pv().lat = [];
% update lat lon
for i = 1:length(pv)
	% check for .1.2.3 and get rid of it from the bus name
	z = pv(i).bus1 == '.'; 
	if ~isempty(find(z,1))
		z = find(z,1,'first');
		z = pv(i).bus1(1:z-1);
	else
		z = pv(i).bus1;
	end
	[x, y] = ismember(z,bl.id);
	pv(i).lon_psudo = bl.coord(y,1);
	pv(i).lat_psudo = bl.coord(y,2);
	pv(i).lon = polyval2(p2lon,pv(i).lon_psudo,pv(i).lat_psudo);
	pv(i).lat = polyval2(p2lat,pv(i).lon_psudo,pv(i).lat_psudo);
end

%% add pv to the footprint
dp = siDeployment('PointLoma_340pvs');
dy = size(dp.ground.longitude,1);
dx = size(dp.ground.longitude,2);

% just take ground map and add PV system to it
gm = dp;
limx = [min(gm.ground.longitude(1,:)) max(gm.ground.longitude(1,:))];
limy = [min(gm.ground.latitude(:,1)) max(gm.ground.latitude(:,1))];

% Locate pixel location for PV system
pv().pixelx = [];
pv().pixely = [];
% dummy index map for calculating the radius
gmidx = repmat(1:dx,dy,1);
gmidy = repmat([1:dy]',1,dx);
% generate the footprint based on log PV size
% reset footprint
gm.footprint.GHI = nan(dy,dx);

%% use typical W/m^2 rating - real scale with rating 100W/m^2 panels and 95% density
% reset footprint
gm.footprint.GHI = nan(dy,dx);
cmap = colormap(jet(256));
pv().dlat = [];
pv().dlon = [];
pv().area = [];
% lon and lat map
l = [limy(1): (limy(2)-limy(1))/(dy-1) : limy(2)]';
latmap = repmat(l,1,dx);
l = [limx(1): (limx(2)-limx(1))/(dx-1) : limx(2)];
lonmap = repmat(l,dy,1);
%
for i = 1:length(pv)
	% PV location
	pv(i).pixelx = round( dx * (pv(i).lon - limx(1)) / (limx(2)-limx(1)) );
	pv(i).pixely = round( dy * (pv(i).lat - limy(1)) / (limy(2)-limy(1)) );
	
	% get pv sizing in real coordinates
	pv_ = pvsystem;
    pv_.cenlat = pv(i).lat;
    pv_.cenlon = pv(i).lon;
    pv_.kVA = pv(i).kVA;
    pv_.density = 0.95;
    pv_.outputrating = 100;
    pv_.calArea;
    
	pv(i).dlat = pv_.dlat;
	pv(i).dlon = pv_.dlon;
	pv(i).area = pv_.area;
	% find the area that in the dlat, dlon specified and convert to pixels
	x =  ( abs(latmap-pv(i).lat) < pv(i).dlat/2 ) & ( abs(lonmap-pv(i).lon) < pv(i).dlon/2 );
	if any(x(:))
		gm.footprint.GHI(x) = i;
	else
		gm.footprint.GHI(pv(i).pixely,pv(i).pixelx) = i;
	end
	
% 	gm.DEMROES_X(i) = pv(i).pixelx;
% 	gm.DEMROES_Y(i) = pv(i).pixely;
end
gm.pvsystem = pv;
figure,imagesc(gm.footprint.GHI);
set(gca,'fontsize',20);
xlabel('X, [pixel]','fontsize',20);
ylabel('Y, [pixel]','fontsize',20);
grid off; box on;set(gcf,'color','w')
set(gcf,'Position',[52 5 1309 940]);
set(gca,'ydir','normal');
% GHI names
for i = 1:length(pv)
	gm.footprint.GHInames{i} = pv(i).Name;
end

%% remove the sensor field if not needed anymore
if isfield(gm.footprint,'sensor')
    gm.footprint = rmfield(gm.footprint,'sensornames');
    gm.footprint = rmfield(gm.footprint,'sensor');
    gm.footprint = rmfield(gm.footprint,'sensorpoints');
end
%% change name fields from GHI to PV 
gm.footprint.PVnames = gm.footprint.GHInames;
gm.footprint.PV = gm.footprint.GHI;
gm.data_type = {'PV'};
gm.footprint = rmfield(gm.footprint,'GHInames');
gm.footprint = rmfield(gm.footprint,'GHI');
%% save back to original deployment folder
deploymentdir = gm.source; ground = gm.ground; design = gm.design; footprint = gm.footprint;
footprint.pvsystem = gm.pvsystem;
save([deploymentdir '/ground.mat'], '-struct', 'ground');
save([deploymentdir '/design.mat'], '-struct', 'design');
save([deploymentdir '/footprint.mat'], '-struct', 'footprint');
% save other parameters in the deployment.conf file
target.data_type = gm.data_type;
writeConf(target,[deploymentdir '/deployment.conf'],1);
msgbox(sprintf('Successfully written to file %s', deploymentdir));
