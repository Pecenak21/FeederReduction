function [result_dir] = plotSimComp( r1, r2, result_dir)
%DAILYPLOT  Function to plot results of dssSimulation
% INPUT:
%   - dataPV: monitors obtained with dssSimulation with pv
%   - dataNoPV: monitors obtained with dssSimulation with pv
%   - result_dir: Path to where to save figures
%   - SimName: simulation name
%   - plotVoltageProfile:   * 0 = do not plot Voltage profile
%                           * 1 = plot and save Voltage profile as png and
%                           fig
%   - plotTime: (optionnal) time in hours at which you want to plot the voltage profil.
%               can be an array with different time.

conf = getConf;
if ~exist('result_dir','var') || isempty(result_dir)
    result_dir = [conf.outputDir '/plot'];
end
if ~exist(result_dir,'dir'), mkdir(result_dir); end

steps = length(r2.VoltMaxMin);
xrange=0:24/steps:24-24/steps;
disp('Plotting results');

%% Max/Min Voltage
f = figure; hold on;
h = plot(xrange, r1.VoltMaxMin(:,1),xrange, r2.VoltMaxMin(:,1),...
    xrange, r1.VoltMaxMin(:,2),xrange, r2.VoltMaxMin(:,2));
legend(h,{'Max (PV)', 'Max (No PV)','Min (PV)', 'Min (No PV)'},'Location','Best');
title('Maximum and Minimum Bus Voltages')
xlabel('Time of Day')
ylabel('Voltages (pu)')
grid on;box on;
xlim([min(xrange) max(xrange)])
saveas(f,[result_dir '/' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' get(get(gca,'title'),'string'), '.fig'], 'fig')

%% Total Loss
% f = figure('units','normalized','outerposition',[0 0 1 1])('visible','off'); hold on;
f = figure; hold on;
h = plot(xrange, r1.TotalLoss(:,1)/1000, xrange, r2.TotalLoss(:,1)/1000,...
    xrange, r1.TotalLoss(:,2)/1000, xrange, r2.TotalLoss(:,2)/1000);
legend(h,{'kW (PV)', 'kW (No PV)', 'kVar (PV)', 'kVar (No PV)'},'Location','Best')
title('Total Losses')
xlabel('Time of Day')
xlim([min(xrange) max(xrange)])
ylabel('Losses')
grid on;box on;
saveas(f,[result_dir '/' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' get(get(gca,'title'),'string'), '.fig'], 'fig')

%% Line Loss
% f = figure('visible','off'); hold on;
f = figure; hold on;
h = plot(xrange, r1.LineLoss(:,1)/1000,  xrange, r2.LineLoss(:,1)/1000,...
    xrange, r1.LineLoss(:,2)/1000,  xrange, r2.LineLoss(:,2)/1000);
legend(h,{'kW (PV)', 'kW (No PV)', 'kVar (PV)', 'kVar (No PV)'},'Location','Best')
title('Total Line Losses')
xlabel('Time of Day')
xlim([min(xrange) max(xrange)])
ylabel('Losses')
grid on;box on;
saveas(f,[result_dir '/' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' get(get(gca,'title'),'string'), '.fig'], 'fig');

%% tap operations
if isfield(r1,'tapPos')
    f1=figure(101); f2=figure(102);
    for j = 1:length(r1.tapPos.transformer)
        fn1 = sprintf('%s/%s%03d_%s',result_dir,'Tranx',j,'TapPos');
        fn2 = sprintf('%s/%s%03d_%s',result_dir,'Tranx',j,'Volt');
%         if ~exist([fn1 '.png'],'file') || ~exist([fn1 '.fig'],'file') || ~exist([fn2 '.png'],'file') || ~exist([fn2 '.fig'],'file')
            % tranx tap position
            figure(101); clf; h1 = plot(r1.time,r1.tapPos.transformer(j).pos,r2.time,r2.tapPos.transformer(j).pos);datetick; grid on;box on; hold off;
            legend(h1,{'S1-Winding 1','S1-Winding 2','S1-Winding 3','S2-Winding 1','S2-Winding 2','S2-Winding 3'},'Location','Best');
            xlabel('Time, [HH:MM]'), ylabel(['Transformer ' num2str(j) '''s Tap Position, [-]']);
            
            % save
            saveFigure(f1,[fn1 '.png']); saveFigure(f1,[fn1 '.fig']);
            
            % tranx voltage
            figure(102); clf; h2 = plot(r1.time,r1.tapPos.transformer(j).V,r2.time,r2.tapPos.transformer(j).V);datetick; grid on;box on; hold off;
            legend(h2,{'S1-Winding 1','S1-Winding 2','S1-Winding 3','S2-Winding 1','S2-Winding 2','S2-Winding 3'},'Location','Best')
            xlabel('Time, [HH:MM]'), ylabel(['Transformer ' num2str(j) '''s Voltage, [p.u.]']);
            % save
            saveFigure(f2,[fn2 '.png']); saveFigure(f2,[fn2 '.fig']);
%         end
    end
end

%% Total Power
f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
h = plot(xrange, r1.TotalPower(:,1),  xrange, r2.TotalPower(:,1),...
    xrange, r1.TotalPower(:,2),  xrange, r2.TotalPower(:,2));
legend(h,{'MW (PV)', 'MW (No PV)', 'MVar (PV)', 'MVar (No PV)'},'Location','Best')
title('Total Power')
xlabel('Time of Day')
xlim([min(xrange) max(xrange)])
ylabel('Power')
grid on;box on;
saveas(f,[result_dir '/' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' get(get(gca,'title'),'string'), '.fig'], 'fig')

%% PV production
f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
h = plot(xrange, -r1.TotalPower(:,1)+r2.TotalPower(:,1),...
    xrange, -r1.TotalPower(:,2)+ r2.TotalPower(:,2));
legend(h,{'MW PV prod', 'difference in MVar'},'Location','Best')
title('PV prod')
xlabel('Time of Day')
xlim([min(xrange) max(xrange)])
ylabel('Power')
grid on;box on;
saveas(f,[result_dir '/' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' get(get(gca,'title'),'string'), '.fig'], 'fig')

end
