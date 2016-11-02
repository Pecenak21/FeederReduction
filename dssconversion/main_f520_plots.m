%% Script to make some plots for report
%% Plot for power output on a clear day
% load clear day data
load('data/20121214_clrSkyGHI.mat');
% load data on clear day (Scenario 2, case 6 for example)
load('scenario2_case6.mat');
% use systems with rated kVA output of 3, 6 and 33kW to plot for
%% demonstration
figure, hold on
set(gcf,'position',[50 50 1000 800])
colors = {'g','r','b','m','b','k'};
x = 0;

for i = [length(c.pvsystem)-10, length(c.pvsystem)-14, length(c.pvsystem)]
	x = x +1;
	plot( ghi.time , ghi.ghi/10^3.*c.pvsystem(i).kVA , colors{x} ,'linewidth',2);
end
datetick('x','HH:MM')
set(gca,'fontsize',20)
ylabel('Power, [kW]','fontsize',20);
xlabel('Time (PST) [HH:MM]','fontsize',20);
xlim([datenum([2012 12 14 06 00 00]) datenum([2012 12 14 18 00 00])]); 
% ylim([0 1000]);
legend('Site 1 (3kW)','Site 15 (6kW)','Site 45 (33kW)');

set(gcf,'color','w')
box on, grid on

%% Fix clear sky curve
figure, set(gcf,'position',[50 50 1000 800]);
plot( ghi.ghi);
y = ghi.ghi; y (1:850) = nan; y(1967:end) = nan;
plot(y);
x = ghi.time - min(ghi.time(:));
[fitresult, gof] = createFit(x,y);
%% extrapolate
y2 = fitresult(x);y2(1:500) = 0; y2(2500:end)=0; y2(y2<0) = 0;
figure, plot(ghi.time,y2,ghi.time,ghi.ghi);

%% fix ghi in clear sky data
ghi.ghi = y2;
save('data/20121214_clrSkyGHI.mat','ghi','pos');

%% plot aggregated load + PV output -> net load

fn = dir('sce*.mat');
for j = 1:length(fn)
    %%
    load(fn(j).name);
    o = [];
    for i = 1:length(c.pvsystem)
        cn = c.pvsystem(i).daily;
        [id, id] = ismember(cn,{c.loadshape.Name}');
        if i == 1
            o = c.loadshape(id).Mult*c.pvsystem(i).kVA;
        else
            o = o + c.loadshape(id).Mult*c.pvsystem(i).kVA;
        end
    end

    % load
    l = [];
    for i = 1:length(c.load)
        cn = c.load(i).Daily;
        [id, id] = ismember(cn,{c.loadshape.Name}');
        if i == 1
            l = c.loadshape(id).Mult*c.load(i).Kw;
        else
            l = l + c.loadshape(id).Mult*c.load(i).Kw;
        end
    end

    %%
    figure(1), set(gcf,'position',[100,100,1100,800]); set(gcf,'color','w');
    if length(l) ==2880
        dt = 30/3600/24;
    elseif length(l) == 24
        dt = 1/24;
    else
        error;
    end
    t = datenum([2012 12 14])+dt:dt:datenum([2012 12 15 00 00 00]);
    plot(t,l/1000,'b','linewidth',2); hold on; plot(t,o/1000,'r','linewidth',2)
    plot(t,(l-o)/1000,'g','linewidth',2); hold off
    datetick
    xlabel('PST Time, HH:MM','fontsize',20);
    xlim([datenum([2012 12 14 00 00 00]) datenum([2012 12 15 00 00 00])]); 
    ylabel('Power, MW','fontsize',20);
    set(gca,'fontsize',20);
    legend('Aggregated Load','Aggregated PV Output','Net Load','location','NorthWest');
    box on; grid on;
    p = 'figure_f520/aggLoad_PV';
    if ~exist(p,'dir')
        mkdir(p);
    end
    %
%     cdata = getframe(gcf); cdata = cdata.cdata;
%     imwrite(cdata,[p '/aggLoadPV_' fn(j).name(1:end-4) '.png']);
    %%
    %saveas(gcf,[p '/aggLoadPV_' fn(j).name(1:end-4) '.png']);
    saveas(gcf,[p '/aggLoadPV_' fn(j).name(1:end-4) '.fig']);
end
