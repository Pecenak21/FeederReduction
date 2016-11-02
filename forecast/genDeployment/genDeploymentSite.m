%% Generate deployment configuration with ground map information for Fallbrook feeder 520
% load feeder520 circuit data. This data is generated from convert2dss tool
c = load('data/feeder520.mat','c'); c = c.c;
% we're only interested in bus location
bl = c.buslist;

% load from disk
load('SDGE_feeders.mat');
figure, hold on
f = f520;
for i = 1:length(f)
	plot([f(i).X]',[f(i).Y]');
end

% let's manually match some buses from psudo-coords to real geographic coords
% do this using circuitvisualizer from convert2dss package and plotting
% function in SDGE_feeders.m
% bus/node: 05201643A (2MW pv system) -> real coords: (x,y) or (lon,lat): [-117.3381 33.4492]
[x, y] = ismember('05201643',bl.id);
i = 1;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-117.3381 33.4492];

% bus/node: 05202234 (Capcitor.520121CW) -> real coords: (x,y) or (lon,lat): [-117.2507 33.4469]
[x, y] = ismember('05202234',bl.id);
i = 2;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-117.2507 33.4469];

% bus/node: 0520 (circuit/substation) -> real coords: (x,y) or (lon,lat): [-117.2368 33.3846]
[x, y] = ismember('0520',bl.id);
i = 3;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-117.2368 33.3846];

% bus/node: 05201314 (capacitor.520_1228CF) -> real coords: (x,y) or (lon,lat): [-117.3019 33.4178]
[x, y] = ismember('05201314',bl.id);
i = 4;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord = [-117.3019 33.4178];

% bus/node: 05202543 (small pv) -> real coords: (x,y) or (lon,lat): [-117.2956 33.4319]
[x, y] = ismember('05202543',bl.id);
i = 5;
bus(i).id = bl.id(y);
bus(i).psudo_coord = bl.coord(y,:);
bus(i).coord =[-117.2956 33.4319];

%% Fit 1D polyfit find x/lon equation. The 2D fitting works much better.
% for i = 1:length(bus)
% 	a(i) = bus(i).psudo_coord(1);
% 	b(i) = bus(i).coord(1);
% end
% figure, plot(a,b,'x-')
% % linear relation
% [plon,Slon] = polyfit(a,b,1);
% 
% % find x/lon equation
% for i = 1:length(bus)
% 	a(i) = bus(i).psudo_coord(2);
% 	b(i) = bus(i).coord(2);
% end
% figure, plot(a,b,'x-')
% % linear relation
% [plat,Slat] = polyfit(a,b,1);
% 
% %% Generate other buses from the list
% buslist(length(bl.coord)) = bus(1);
% for i = 1:length(buslist)
% 	buslist(i).id = bl.id(i);
% 	buslist(i).psudo_coord = bl.coord(i,:);
% 	buslist(i).coord(1) = polyval(plon,bl.coord(i,1),Slon);
% 	buslist(i).coord(2) = polyval(plat,bl.coord(i,2),Slat);
% end
% 
% %% plot again to compare
% load('SDGE_feeders.mat');
% figure, hold on
% f = f520;
% for i = 1:length(f)
% 	plot([f(i).X]',[f(i).Y]');
% end
% %
% hold on
% for i = 1:length(buslist)
% 	plot(buslist(i).coord(1), buslist(i).coord(2),'.r');
% end

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

%% plot again to compare
% load('SDGE_feeders.mat');
% figure, hold on
% f = f520;
% for i = 1:length(f)
% 	plot([f(i).X]',[f(i).Y]');
% end
% %
% % hold on
% % for i = 1:length(buslist)
% % 	plot(buslist(i).coord(1), buslist(i).coord(2),'.r');
% % end
% 
% hold on
% for i = 1:length(buslist2)
% 	plot(buslist2(i).coord(1), buslist2(i).coord(2),'.g');
% end
% set(gca,'ydir','normal');

%% Load PV data. The data load won't be in object types like in convert2dss package but we can deal with that.
c = load('data/feeder520_wPV.mat','c'); c = c.c;
pv = [c.pvsystem.data];

% clean up data
x = setdiff(fields(pv),{'kVA','bus1','Name'});
pv = rmfield(pv,x);
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

%% Plot PV with log scale for power output
load('SDGE_feeders.mat');
figure, hold on
f = f520;
limx = [-117.35 -117.21];
limy = [33.375 33.48];
for i = 1:length(f)
	plot([f(i).X]',[f(i).Y]','k');
end
grid off; box on;set(gcf,'color','w')
set(gcf,'Position',[52 5 1309 940]);
xlim(limx); ylim(limy);

% linear scale
cmap = colormap(jet(256));
for i = 1:length(pv)
		plot([pv(i).lon]',[pv(i).lat]','.r','markersize',log(pv(i).kVA)*15,'linewidth',2,...
			'color',cmap(ceil( log10(pv(i).kVA)/log10(1000)*255 ),:) );
end
xlabel('Longitude','fontsize',20);
ylabel('Latitude','fontsize',20);
set(gca,'fontsize',20);
set(gca,'ydir','normal');
% Rescale color bar
d = log10([pv.kVA]);
mn = min(d(:));
rng = max(d(:))-mn;
d = 1+63*(d-mn)/rng; % Self scale data
hC = colorbar;
L = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000];
% Choose appropriate
% or somehow auto generate colorbar labels
l = 1+255*(log10(L)-mn)/rng; % Tick mark positions
set(hC,'Ytick',l,'YTicklabel',L,'fontsize',20);
ylabel(hC,'PV system''s rated output, [kW]','fontsize',20);

%% PV location on the UCSD grid
% generate ground map from PV map
% figure, hold on
% grid off; box on;set(gcf,'color','w')
% set(gcf,'Position',[52 5 1309 940]);
% xlabel('Longitude','fontsize',20);
% ylabel('Latitude','fontsize',20);
% set(gca,'fontsize',20);
% limx = [-117.35 -117.21];
% limy = [33.375 33.48];
% xlim(limx); ylim(limy); 
% colormap(jet(256));
% hC = colorbar;

% get UCSD footprint
ucsd = siDeployment('UCSD');
dy = size(ucsd.footprint.GHI,1);
dx = size(ucsd.footprint.GHI,2);

% initialize groundmap
l = [limy(1): (limy(2)-limy(1))/(dy-1) : limy(2)]';
gm.ground.latitude = repmat(l,1,dx);
l = [limx(1): (limx(2)-limx(1))/(dx-1) : limx(2)];
gm.ground.longitude = repmat(l,dy,1);

% just take UCSD ground map and add PV system to it
gm = ucsd;
gm.name = 'Fallbrook Feeder';

% Locate pixel location for PV system
% new fields
pv().pixelx = [];
pv().pixely = [];
% dummy index map for calculating the radius
gmidx = repmat(1:dx,dy,1);
gmidy = repmat([1:dy]',1,dx);
% generate the footprint based on log PV size
% reset footprint
gm.footprint.GHI = nan(dy,dx);
cmap = colormap(jet(256));
for i = 1:length(pv)
	% PV location
	pv(i).pixelx = round( dx * (pv(i).lon - limx(1)) / (limx(2)-limx(1)) );
	pv(i).pixely = round( dy * (pv(i).lat - limy(1)) / (limy(2)-limy(1)) );
	
% 	cid = ceil( log10(pv(i).kVA)/log10(1000)*255 );
	% to make sure each PV system has a unique index value
% 	while ismember(cid,gm.footprint.GHI)
% 		if cid < 255
% 			cid = cid + 1;
% 		else
% 			cid = cid - 1;
% 		end
% 	end
% 	pv.idsize = cid;
% 	gm.footprint.GHI( pv(i).pixely, pv(i).pixelx ) = cid;
	% PV size
	s = log(pv(i).kVA)*5;
	x = ( (gmidx-pv(i).pixelx).^2 + (gmidy-pv(i).pixely).^2 ).^.5 < s;
	gm.footprint.GHI(x) = i;
end
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

%% generate the footprint based on typical W/m^2 rating - Linear scaling
% reset footprint
gm.footprint.GHI = nan(dy,dx);
cmap = colormap(jet(256));
for i = 1:length(pv)
	% PV location
	pv(i).pixelx = round( dx * (pv(i).lon - limx(1)) / (limx(2)-limx(1)) );
	pv(i).pixely = round( dy * (pv(i).lat - limy(1)) / (limy(2)-limy(1)) );
	
% 	cid = ceil( log10(pv(i).kVA)/log10(1000)*255 );
	% to make sure each PV system has a unique index value
% 	while ismember(cid,gm.footprint.GHI)
% 		if cid < 255
% 			cid = cid + 1;
% 		else
% 			cid = cid - 1;
% 		end
% 	end
% 	pv.idsize = cid;
% 	gm.footprint.GHI( pv(i).pixely, pv(i).pixelx ) = cid;
	% PV size
	s = sqrt(pv(i).kVA)*2;
	x = ( (gmidx-pv(i).pixelx).^2 + (gmidy-pv(i).pixely).^2 ).^.5 < s;
	gm.footprint.GHI(x) = i;
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
%% save to Fallbrook_PV_linear_scaling
filen = 'Fallbrook/footprint.mat';
fp = load(filen);
% change UCSD footprint to Fallbrook footprint
fp.GHI = flipdim(gm.footprint.GHI,1);
fn = fields(fp)';
%
for i = 1:length(fn);
	if i == 1; save(filen,'-struct','fp',fn{i});
	else
		save(filen,'-struct','fp',fn{i},'-append');
	end
end

%% generate the footprint based on typical W/m^2 rating - real scale with rating 100W/m^2 panels and 95% density
% reset footprint
gm.footprint.GHI = nan(dy,dx);
cmap = colormap(jet(256));
pv().dlat = [];
pv().dlon = [];
pv().area = [];
pv().onesize = [];
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
	pv_ = PVSystem(pv(i).lat,pv(i).lon, pv(i).kVA, .95, 100);
	pv(i).dlat = pv_.dlat;
	pv(i).dlon = pv_.dlon;
	pv(i).area = pv_.area;
	pv(i).onesize = pv_.squaresize;
	% find the area that in the dlat, dlon specified and convert to pixels
	x =  ( abs(latmap-pv(i).lat) < pv(i).dlat/2 ) & ( abs(lonmap-pv(i).lon) < pv(i).dlon/2 );
	if any(x(:))
		gm.footprint.GHI(x) = i;
	else
		gm.footprint.GHI(pv(i).pixely,pv(i).pixelx) = i;
	end
	
	gm.DEMROES_X(i) = pv(i).pixelx;
	gm.DEMROES_Y(i) = pv(i).pixely;
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

%% save to virtual Fallbrook feeder based on UCSD ground map
filen = 'Fallbrook/footprint.mat';
fp = load(filen);
% change UCSD footprint to Fallbrook footprint
fp.GHI = flipdim(gm.footprint.GHI,1);
fn = fields(fp)';
%
for i = 1:length(fn);
	if i == 1; save(filen,'-struct','fp',fn{i});
	else
		save(filen,'-struct','fp',fn{i},'-append');
	end
end

%% map size
map.dxkm = greatCircleDistance(limy(1),limx(1),limy(1),limx(2));
map.dykm = greatCircleDistance(limy(1),limx(1),limy(2),limx(1));

%% package up the data to different PV profiles
load('/home/andu/code/forecast/FallbrookPVProfilesRealScale_USI1_2/forecast_001.mat');

dirname = 'FallbrookPVProfilesRealScale_USI1_2';
fl = dir(fullfile(dirname,'*.mat'));
for i = 1:length(fl)
	load( [dirname '/' fl(i).name ] );
end

pv().gi = [];
pv().time = [];
pv().power = [];

tlength = size(time,1);
forecastIndex = 1; % nowcast
efficiency = .1; % 10% effi or 100W/m^2 for 1000W/m^2 irradiance
for j = 1:length(pv)
	pv(j).time = time(:,1);
	pv(j).gi = nan(tlength,1);
	pv(j).power = nan(tlength,1);
	for i = 1:tlength
		pv(j).gi(i) = station(i,j).gi(forecastIndex);
		pv(j).power(i) = pwrGIToPower(pv(j).gi(i),1000)*pv(j).area*.95*efficiency;
	end
end

save([dirname '/pv.mat'],'pv');

%% plot GI 2 hours with 1MW site
figure, hold on
colors = {'g','r','b','m','b','k'};
x = 0;
for i = [1 15 39]
	x = x +1;
	plot( pv(i).time - 8/24 , pv(i).gi , colors{x} ,'linewidth',2);
end
datetick('x','HH:MM')
set(gca,'fontsize',20)
ylabel('GHI, [W m^{-2}]','fontsize',20);
xlabel('Time (PST) [HH:MM]','fontsize',20);
xlim([datenum([2012 12 14 10 00 00]) datenum([2012 12 14 12 00 00])]); 
ylim([0 1000]);
legend('Site 1 (3kW)','Site 15 (6kW)','Site 39 (1MW)');
set(gcf,'position',[50 50 1000 800])
set(gcf,'color','w')
box on, grid on

%% plot GI whole day with 1MW site
figure, hold on
colors = {'g','r','b','m','b','k'};
x = 0;
for i = [1 15 39]
	x = x +1;
	plot( pv(i).time - 8/24 , pv(i).gi , colors{x} ,'linewidth',2);
end
datetick('x','HH:MM')
set(gca,'fontsize',20)
ylabel('GHI, [W m^{-2}]','fontsize',20);
xlabel('Time (PST) [HH:MM]','fontsize',20);
xlim([datenum([2012 12 14 08 00 00]) datenum([2012 12 14 16 00 00])]); 
ylim([0 1000]);
legend('Site 1 (3kW)','Site 15 (6kW)','Site 39 (1MW)');
set(gcf,'position',[50 50 1000 800])
set(gcf,'color','w')
box on, grid on

%% plot power 2 hours with 1MW site
figure, hold on
colors = {'g','r','b','m','b','k'};
x = 0;
for i = [1 15 39]
	x = x +1;
	plot( pv(i).time - 8/24 , pv(i).power/10^6 , colors{x} ,'linewidth',2);
end
datetick('x','HH:MM')
set(gca,'fontsize',20)
ylabel('Power, [MW]','fontsize',20);
xlabel('Time (PST) [HH:MM]','fontsize',20);
xlim([datenum([2012 12 14 10 00 00]) datenum([2012 12 14 12 00 00])]); 
% ylim([0 1000]);
legend('Site 1 (3kW)','Site 15 (6kW)','Site 39 (1MW)');
set(gcf,'position',[50 50 1000 800])
set(gcf,'color','w')
box on, grid on

%% plot power whole day with 1MW site
figure, hold on
colors = {'g','r','b','m','b','k'};
x = 0;
for i = [1 39 15]
	x = x +1;
	plot( pv(i).time - 8/24 , pv(i).power/10^6 , colors{x} ,'linewidth',2);
end
datetick('x','HH:MM')
set(gca,'fontsize',20)
ylabel('Power, [MW]','fontsize',20);
xlabel('Time (PST) [HH:MM]','fontsize',20);
xlim([datenum([2012 12 14 08 00 00]) datenum([2012 12 14 16 00 00])]); 
% ylim([0 1000]);
legend('Site 1 (3kW)','Site 39 (1MW)','Site 15 (6kW)');
set(gcf,'position',[50 50 1000 800])
set(gcf,'color','w')
box on, grid on

%% plot power 2 hours with 33kw site
figure, hold on
colors = {'g','r','b','m','b','k'};
x = 0;
for i = [1 15 45]
	x = x +1;
	plot( pv(i).time - 8/24 , pv(i).power/10^3 , colors{x} ,'linewidth',2);
end
datetick('x','HH:MM')
set(gca,'fontsize',20)
ylabel('Power, [kW]','fontsize',20);
xlabel('Time (PST) [HH:MM]','fontsize',20);
xlim([datenum([2012 12 14 10 00 00]) datenum([2012 12 14 12 00 00])]); 
% ylim([0 1000]);
legend('Site 1 (3kW)','Site 15 (6kW)','Site 45 (33kW)');
set(gcf,'position',[50 50 1000 800])
set(gcf,'color','w')
box on, grid on

%% plot power whole day with 33kw site
figure, hold on
colors = {'g','r','b','m','b','k'};
x = 0;
for i = [1 15 45]
	x = x +1;
	plot( pv(i).time - 8/24 , pv(i).power/10^3 , colors{x} ,'linewidth',2);
end
datetick('x','HH:MM')
set(gca,'fontsize',20)
ylabel('Power, [kW]','fontsize',20);
xlabel('Time (PST) [HH:MM]','fontsize',20);
xlim([datenum([2012 12 14 08 00 00]) datenum([2012 12 14 16 00 00])]); 
% ylim([0 1000]);
legend('Site 1 (3kW)','Site 15 (6kW)','Site 45 (33kW)');
set(gcf,'position',[50 50 1000 800])
set(gcf,'color','w')
box on, grid on

%% Plot kt to check
figure, plot(time(:,1)-8/24,kt(:,1),time(:,1)-8/24,kt(:,2),time(:,1)-8/24,kt(:,3),'linewidth',2);
%
ylabel('kt','fontsize',20);
datetick('x','HH:MM')
set(gca,'fontsize',20)
xlabel('Time (PST) [HH:MM]','fontsize',20);
xlim([datenum([2012 12 14 08 00 00]) datenum([2012 12 14 16 00 00])]); 
set(gcf,'position',[50 50 1000 800])
set(gcf,'color','w')
box on, grid on
legend('Thick Cloud','Thin Cloud','Clear');

