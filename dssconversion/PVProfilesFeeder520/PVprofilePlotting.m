load pv.mat
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
