clear
close
% load('R:\infobase\smart Inverters\dynamicSimulations\feeder Reduction\reducedCircuits/Alpine_Reducedcircuit.mat')
% load('R:\infobase\smart Inverters\dynamicSimulations\feeder Reduction\reducedCircuits/PointLoma_Reducedcircuit.mat')
% load('R:\infobase\smart Inverters\dynamicSimulations\feeder Reduction\reducedCircuits/Ramona_Reducedcircuit.mat')
% load('R:\infobase\smart Inverters\dynamicSimulations\feeder Reduction\reducedCircuits/ValleyCenter_Reducedcircuit.mat')
% load('c:\users\zactus\gridIntegration\results\Artesian_wpv_existing.mat')
% load('R:\infobase\smart Inverters\dynamicSimulations\feeder Reduction\reducedCircuits/Fallbrook_Reducedcircuit.mat')
% load('c:\users\zactus\gridIntegration\results\Alpine_s1b.mat')
load('R:/infobase/smart Inverters/dynamicSimulations/feeder Reduction/reducedCircuits/fallbrook_Reducedcircuit.mat')
busName=lower(c.buslist.id);
bus_coord=c.buslist.coord;
line = c.line;
trf=c.transformer;
swt=c.switch;

usName = lower(busName);
% flag=zeros(length(line),2);
figure;
for i = 1:length(line);
	flag(i,1) = find(ismember(busName,lower(regexp(line(i).bus1,'\.','split','once'))));
	flag(i,2) = find(ismember(busName,lower(regexp(line(i).bus2,'\.','split','once'))));
	x_tmp(1,1)= bus_coord(flag(i,1),1);  x_tmp(1,2)= bus_coord(flag(i,2),1);
	y_tmp(1,1) = bus_coord(flag(i,1),2); y_tmp(1,2) =bus_coord(flag(i,2),2);
	h(1)=plot(x_tmp,y_tmp,'Color',[.65 .65 .65],'linewidth',3);
	hold on;
% 	  text(x_tmp(1,1),y_tmp(1,1),busName(flag(i,1)))
% 		text(x_tmp(1,2),y_tmp(1,2),busName(flag(i,2)))
end

% for i=1:length(swt)
% 	LinInd=find(ismember([line(:).Name],regexprep(swt(i).SwitchedObj,'Line.','')));
% 	flag(i,1) = find(ismember(busName,lower(regexp(line(LinInd).bus1,'\.','split','once'))));
% 	flag(i,2) = find(ismember(busName,lower(regexp(line(LinInd).bus2,'\.','split','once'))));
% 	x_tmp(1,1)= bus_coord(flag(i,1),1);  x_tmp(1,2)= bus_coord(flag(i,2),1);
% 	y_tmp(1,1) = bus_coord(flag(i,1),2); y_tmp(1,2) =bus_coord(flag(i,2),2);
% 	plot(x_tmp(1,1),y_tmp(1,1),'*g',x_tmp(1,2),y_tmp(1,2),'*r','linewidth',3);
% 	 text(x_tmp(1,1),y_tmp(1,1),swt(i).Action)
% end
busName=lower(circuit.buslist.id);
bus_coord=circuit.buslist.coord;
line = circuit.line;
trf=circuit.transformer;

for i = 1:length(line);
	flag(i,1) = find(ismember(busName,lower(regexp(line(i).bus1,'\.','split','once'))));
	flag(i,2) = find(ismember(busName,lower(regexp(line(i).bus2,'\.','split','once'))));
	x_tmp(1,1)= bus_coord(flag(i,1),1);  x_tmp(1,2)= bus_coord(flag(i,2),1);
	y_tmp(1,1) = bus_coord(flag(i,1),2); y_tmp(1,2) =bus_coord(flag(i,2),2);
	h(2)=plot(x_tmp,y_tmp,'Color','k','linewidth',3);
	hold on;
	  text(x_tmp(1,1),y_tmp(1,1),busName(flag(i,1)))
		text(x_tmp(1,2),y_tmp(1,2),busName(flag(i,2)))
		text((x_tmp(1,2)+x_tmp(1,1))/2+100,(y_tmp(1,2)+y_tmp(1,1))/2+100,num2str(i))
end

% T=[39  325 382 462 555 549 293 319 100 614];
% h(2)=plot(bus_coord(T,1),bus_coord(T,2),'gx','linewidth',1.5,'MarkerFaceColor','g','markersize',14,'markeredgecolor','g');
% 
% if isfield(circuit,'pvsystem')
% 	pv=circuit.pvsystem;
% 	for i=1:length(pv)
% 		bus1=regexp(pv(i).bus1,'\.','split','once');
% 		flag(i,1) = find(ismember(lower(busName),lower(bus1{1})));
% 		h(3)=plot(bus_coord(flag(i,1),1),bus_coord(flag(i,1),2),'g>','linewidth',1.5,'MarkerFaceColor','g','markersize',4,'markeredgecolor','g');
% 	end
% end

%Plot USI

% % hold on;plot(usi.position.longitude,usi.position.latitude,'rx')
% % r= 500/110899;
% % %plot radius
% % th = 0:pi/50:2*pi;
% % xunit = r * cos(th) + usi.position.longitude;
% % yunit = r * sin(th) + usi.position.latitude;
% % plot(xunit, yunit);
% % plot_google_map('MapType', 'hybrid');

% if isfield(circuit,'capacitor')
% 	cap=circuit.capacitor;
% 	for i=1:length(cap)
%     flag(i,1) = find(ismember(busName,lower(regexp(cap(i).bus1,'\.','split','once'))));
% 	h(5)=plot(bus_coord(flag(i,1),1),bus_coord(flag(i,1),2),'cs','linewidth',1.5,'MarkerFaceColor','c','markersize',8,'markeredgecolor','c');
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
% load('c:\users\zactus\gridIntegration\results\Alpine_reduced.mat')
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
% legend(h([1 2 3 4 5]),'Feeder Lines','Reduced Lines','pv systems','Transformer','Capacitor')
set(gca,'visible','off')
set(gca,'XtickLabel',[],'YtickLabel',[]);