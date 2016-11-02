function [o, data] = plotFeederProfile(c)

if ~isa(c,'COM.OpendssEngine_dss')
	o = dssget(c,[],[],[],1);
else
	o = c;
end
	
%% Plot voltage 
figure,plotVoltageProfile(c);

%% losses

%% plot power profile
bnm = o.ActiveCircuit.AllBusNames; %bus name
bdist = o.ActiveCircuit.AllBusDistances'; 
pw = dssgetval(o,c.line,'powers'); 
pwbus2 = nansum(pw(:,6+[1,3,5]),2) + 1i*nansum(pw(:,6+[2,4,6]),2);
pwbus1 = nansum(pw(:,[1,3,5]),2) + 1i*nansum(pw(:,[2,4,6]),2);
[y1, bus1id] = ismember(bnm,cleanBus({c.line.bus1})); bus1id = bus1id(bus1id>0);
[y2, bus2id] = ismember(bnm,cleanBus({c.line.bus2})); bus2id = bus2id(bus2id>0);
p = nan(length(bnm),1);
p(y2) = pwbus2(bus2id);
p(y1) = pwbus1(bus1id);
figure, plot(  bdist, abs(p),'.', bdist, real(p), '.', bdist, imag(p), '.' );
legend({'Total Apparent Power [kVA]','Active Power [kW]','Reactive Power [kVar]'});
grid on; box on; set(gcf,'color','w'); %ylim([.9 1.1]);
set(gca,'fontsize',20);
xlabel('Distance'); 
ylabel('Power [kX]'); 
title('Power Profile','fontsize',20); 
% set(gca,'xtick',[0:2:20])

%% plot current profile
cu = dssgetval(o,c.line,'currents'); 
cbus2 = nansum(cu(:,6+[1,3,5]),2) + 1i*nansum(cu(:,6+[2,4,6]),2);
cbus1 = nansum(cu(:,[1,3,5]),2) + 1i*nansum(cu(:,[2,4,6]),2);
cc = nan(length(bnm),1);
cc(y2) = cbus2(bus2id);
cc(y1) = cbus1(bus1id);
figure, plot(  bdist, abs(cc),'.');
%legend({'Total Apparent Power [kVA]','Active Power [kW]','Reactive Power [kVar]'});
grid on; box on; set(gcf,'color','w'); %ylim([.9 1.1]);
set(gca,'fontsize',20);
title('Current Profile','fontsize',20);
xlabel('Distance'); 
ylabel('Current [A]'); 
% set(gca,'xtick',[0:2:20])

%% output data
data.busName = bnm;
data.busDist = bdist;
data.busPower = p;
data.busCurrent = cc;
end