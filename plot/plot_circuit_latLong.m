clear
close
% load('c:\users\zactus\gridIntegration\results\PointLoma_wpv_existing.mat')
load('c:\users\zactus\gridIntegration\results\ValleyCenter_wpv_existing.mat')
busName=c.buslist.id;
bus_coord=c.buslist.coord;
line = c.line;
trf=c.transformer;

if strcmp(lower(c.circuit.Name),'cabrillo')
	%% get bus location in both psudo and lat-lon coordinates
	% psudo buslist location from SDGE data
	% load('SDGE_feeders.mat');
	% figure, hold on
	% f = f480;
	% for i = 1:length(f)
	% 	plot([f(i).X]',[f(i).Y]');
	% end
	
	% manually map psudo and real coordinates of a few buses
	[x, y] = ismember('04809924',busName);
	i = 1;
	bus(i).id = busName(y);
	bus(i).psudo_coord = bus_coord(y,:);
	bus(i).coord = [-117.2446833, 32.7421833];
	
	% bus/node: 04804780 -> real coords: (x,y) or (lon,lat): [ -117.2564472°, 32.7261500°]
	[x, y] = ismember('04804780',busName);
	i = 2;
	bus(i).id = busName(y);
	bus(i).psudo_coord = bus_coord(y,:);
	bus(i).coord = [-117.2564472, 32.7261500];
	
	% bus/node: 04801426 -> real coords: (x,y) or (lon,lat): [-117.2495222°,32.7082167°]
	[x, y] = ismember('04801426',busName);
	i = 3;
	bus(i).id = busName(y);
	bus(i).psudo_coord = bus_coord(y,:);
	bus(i).coord = [-117.2495222,32.7082167];
	
	% bus/node: 04803609  -> real coords: (x,y) or (lon,lat): [-117.2409222°, 32.7250861°]
	[x, y] = ismember('04803609',busName);
	i = 4;
	bus(i).id = busName(y);
	bus(i).psudo_coord = bus_coord(y,:);
	bus(i).coord = [-117.2409222, 32.7250861];
	
	% bus/node: 04808232 -> real coords: (x,y) or (lon,lat): [-117.2493417°, 32.7303222°]
	[x, y] = ismember('04808232',busName);
	i = 5;
	bus(i).id = busName(y);
	bus(i).psudo_coord = bus_coord(y,:);
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
	buslist2(length(bus_coord)) = bus(1);
	for i = 1:length(buslist2)
		buslist2(i).id = busName(i);
		buslist2(i).psudo_coord = bus_coord(i,:);
		buslist2(i).coord(1) = polyval2(p2lon,bus_coord(i,1),bus_coord(i,2));
		buslist2(i).coord(2) = polyval2(p2lat,bus_coord(i,1),bus_coord(i,2));
	end
	
	usi=siImager('usi_1_9');
elseif strcmp(lower(c.circuit.Name),'valley_center')
	% let's manually match some buses from psudo-coords to real geographic coords
	% do this using circuitvisualizer from convert2dss package and plotting
	% function in SDGE_feeders.m
	% bus/node: 09099435  -> real coords: (x,y) or (lon,lat): [-116.93067500000001°, 33.26428888888889°]
	[x, y] = ismember('09099435',busName);
	i = 1;
	bus(i).id = busName(y);
	bus(i).psudo_coord = bus_coord(y,:);
	bus(i).coord = [-116.93067500000001, 33.26428888888889];
	
	% bus/node: 0909 -> real coords: (x,y) or (lon,lat): [ -117.01785833333334°, 33.23080277777778°]
	[x, y] = ismember('0909',busName);
	i = 2;
	bus(i).id = busName(y);
	bus(i).psudo_coord = bus_coord(y,:);
	bus(i).coord = [ -117.01785833333334, 33.23080277777778];
	
	% bus/node: 09093928 -> real coords: (x,y) or (lon,lat): [-116.98524166666667°,33.220394444444445°]
	[x, y] = ismember('09093928',busName);
	i = 3;
	bus(i).id = busName(y);
	bus(i).psudo_coord = bus_coord(y,:);
	bus(i).coord = [-116.98524166666667,33.220394444444445];
	
	% bus/node: 0909103  -> real coords: (x,y) or (lon,lat): [-116.95702222222222°, 33.268880555555555°]
	[x, y] = ismember('0909103',busName);
	i = 4;
	bus(i).id = busName(y);
	bus(i).psudo_coord = bus_coord(y,:);
	bus(i).coord = [-116.95702222222222, 33.268880555555555];
	
	% bus/node: 09097319 -> real coords: (x,y) or (lon,lat): [-116.95795555555556°, 33.23663055555556°]
	[x, y] = ismember('09097319',busName);
	i = 5;
	bus(i).id = busName(y);
	bus(i).psudo_coord = bus_coord(y,:);
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
	buslist2_909(length(bus_coord)) = bus(1);
	for i = 1:length(buslist2_909)
		buslist2_909(i).id = busName(i);
		buslist2_909(i).psudo_coord = bus_coord(i,:);
		buslist2_909(i).coord(1) = polyval2(p2lon,bus_coord(i,1),bus_coord(i,2));
		buslist2_909(i).coord(2) = polyval2(p2lat,bus_coord(i,1),bus_coord(i,2));
	end
	
	usi=siImager('usi_1_2');
end

busName = lower(busName);
% flag=zeros(length(line),2);
figure;
for i = 1:length(line);
	flag(i,1) = find(ismember(busName,lower(regexp(line(i).bus1,'\.','split','once'))));
	flag(i,2) = find(ismember(busName,lower(regexp(line(i).bus2,'\.','split','once'))));
	x_tmp(1,1)= buslist2(flag(i,1)).coord(1);  x_tmp(1,2)= buslist2(flag(i,2)).coord(1);
	y_tmp(1,1) = buslist2(flag(i,1)).coord(2); y_tmp(1,2) =buslist2(flag(i,2)).coord(2);
	h(1)=plot(x_tmp,y_tmp,'Color',[.65 .65 .65],'linewidth',3);
	hold on;
	% %   text(x_tmp(1,1),y_tmp(1,1),busName(flag(i,1)))
	% % 	text(x_tmp(1,2),y_tmp(1,2),busName(flag(i,2)))
end

% T=[39  325 382 462 555 549 293 319 100 614];
% h(2)=plot(bus_coord(T,1),bus_coord(T,2),'gx','linewidth',1.5,'MarkerFaceColor','g','markersize',14,'markeredgecolor','g');

if isfield(c,'pvsystem')
	pv=c.pvsystem;
	for iiii=1:length(pv)
		bus1=regexp(pv(iiii).bus1,'\.','split','once');
		flag(i+iiii,1) = find(ismember(lower(busName),lower(bus1{1})));
		h(8)=plot(buslist2(flag(i+iiii)).coord(1),buslist2(flag(i+iiii)).coord(2),'g>','linewidth',1.5,'MarkerFaceColor','g','markersize',4,'markeredgecolor','g');
	end
end

%Plot USI

hold on;plot(usi.position.longitude,usi.position.latitude,'rx')
r= 500/110899;
%plot radius
th = 0:pi/50:2*pi;
xunit = r * cos(th) + usi.position.longitude;
yunit = r * sin(th) + usi.position.latitude;
plot(xunit, yunit);
plot_google_map('MapType', 'hybrid');

% for ii=1:length(trf)
% 		buses{ii}=trf(ii).Buses;
% 		Bus1name=regexp(buses{ii}(1),'\.','split');
% 		bus1{ii}=char(Bus1name{1}(1));
% 		Bus2name=regexp(buses{ii}(2),'\.','split');
% 		bus2{ii}=char(Bus2name{1}(1));
% 		flag(i+ii,1) = find(ismember(lower(busName),lower(bus1{ii})));
% 	    flag(i+ii,2) = find(ismember(lower(busName),lower(bus2{ii})));
% 		x_tmp(1,1)= bus_coord(flag(i+ii,1),1);  x_tmp(1,2)= bus_coord(flag(i+ii,2),1);
% 		y_tmp(1,1) = bus_coord(flag(i+ii,1),2); y_tmp(1,2) = bus_coord(flag(i+ii,2),2);
% 		h(3)=plot(x_tmp,y_tmp,'ro','linewidth',1.5,'MarkerFaceColor','b','markersize',8,'markeredgecolor','k');
% % 		text(x_tmp(1,1),y_tmp(1,1),num2str(flagy(ii,1)))
% % 		text(x_tmp(1,2),y_tmp(1,2),num2str(flagy(ii,2)))
% end
%
% if isfield(c,'capacitor')
% 	cap=c.capacitor;
% 	for iii=1:length(cap)
%     flag(i+ii+iii,1) = find(ismember(busName,lower(regexp(cap(iii).bus1,'\.','split','once'))));
% 	h(4)=plot(bus_coord(flag(i+ii+iii,1),1),bus_coord(flag(i+ii+iii,1),2),'cs','linewidth',1.5,'MarkerFaceColor','c','markersize',8,'markeredgecolor','c');
% 	end
% end
%
%
%
% %topo nodes
% topo=[10;27;37;41;161;248;249;312;359;461;545;559];
% 	h(5)=plot(bus_coord(topo,1),bus_coord(topo,2),'rx','linewidth',1.5,'MarkerFaceColor','r','markersize',14,'markeredgecolor','r');
%
% hold on;
% h(6)=plot(bus_coord(end,1),bus_coord(end,2),'bs','markersize',12);
% ylim([236000 247000])
% axis off
% box off
% % h_legend=legend(h([1 2 3 4 5 6]),'Feeder Lines','Choosen Critical Nodes','Transformer','Capacitor','Connecting Nodes','Substation','location','sw');
% % set(h_legend,'FontSize',14);
% % legend(h([1 10000 20000 30000 40000 50000]),'Feeder Lines','Choosen Critical Bus','Transformer','Capacitor','Connecting Buses','Substation','location','sw')
%
% % legend(h([1 10000 20000 30000 40000 50000]),'Feeder Lines','Choosen Critical Bus','Transformer','Capacitor','Connecting Buses','Substation','location','bestoutside')
%
%
% % legend(h([1 8 3aa 4 6]),'Feeder Lines','pvsystems','Transformer','Capacitor', 'substation')
% load('c:\Users\zactus\gridIntegration\results\Alpine_reduced.mat')
% busName=c.buslist.id;
% bus_coord=c.buslist.coord;
% line = c.line;
% trf=c.transformer;
% busName = lower(busName);
%
% for i = 1:length(line);
%     flag(i,1) = find(ismember(busName,lower(regexp(line(i).bus1,'\.','split','once'))));
%     flag(i,2) = find(ismember(busName,lower(regexp(line(i).bus2,'\.','split','once'))));
%     x_tmp(1,1)= bus_coord(flag(i,1),1);  x_tmp(1,2)= bus_coord(flag(i,2),1);
%     y_tmp(1,1) = bus_coord(flag(i,1),2); y_tmp(1,2) = bus_coord(flag(i,2),2);
%     h(7)=plot(x_tmp,y_tmp,'r','linewidth',1.5);
%     hold on;
% end
% h_legend=legend(h([1 2 3 4 5 6 7]),'Feeder Lines','Choosen Critical Bus','Transformer','Capacitor','Connecting Buses','Substation','Reduced Lines')
% set(h_legend,'FontSize',14);