%% plot aggregated load + PV output -> net load for scenario 4

% Vadim: I assume you saved the .mat circuit output files for scenario 4 in
% current directory
d = pwd;
fn = dir([d '/scenario4_*.mat']);
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
    saveas(gcf,[p '/aggLoadPV_' fn(j).name(1:end-4) '.png']);
    saveas(gcf,[p '/aggLoadPV_' fn(j).name(1:end-4) '.fig']);
end