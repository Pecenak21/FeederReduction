function DailyPlot( sim1, sim2, result_dir, SimName, plotVoltageProfile, plotTime )
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
    if ~exist(result_dir,'dir'), mkdir(result_dir); end
end

steps = length(sim2.VoltMaxMin);
xrange=0:24/steps:24-24/steps;
colorCurve = [{'r*-'} , {'g*-'},{'b*-'},{'k*-'},{'y*-'},{'m*-'},{'rx:'},{'gx:'},{'bx:'},{'kx:'}];
ColorLine = [{'k'} , {'b'},{'c'},{'y'},{'m'},{'r'},{'g'}];
disp(['Plotting results for ' SimName]);

%% Fill variables:
% PV:
Volt_MaxMin_PV = sim1.VoltMaxMin;
LossTotal_PV = sim1.TotalLoss;
LossLine_PV = sim1.LineLoss;
TotalPower_PV = sim1.TotalPower;
Event_PV = sim1.Events;
Cap_PV = sim1.Capcontrol;
Reg_PV = sim1.Regulation;
tap2plotPV = sim1.Tap2Plot;
cap2plotPV = sim1.Cap2Plot;
Volt_PV = sim1.Voltage;
DistPV = sim1.Dist;
NodeNamePV=sim1.nodeName;

% No PV:
Volt_MaxMin_NoPV = sim2.VoltMaxMin;
LossTotal_NoPV = sim2.TotalLoss;
LossLine_NoPV = sim2.LineLoss;
Event_NoPV =sim2.Events;
Cap_NoPV = sim2.Capcontrol;
Reg_NoPV = sim2.Regulation;
Volt_NoPV = sim2.Voltage;
TotalPower_NoPV = sim2.TotalPower;
tap2plotNoPV = sim2.Tap2Plot;
cap2plotNoPV = sim2.Cap2Plot;
DistNoPV = sim2.Dist;
NodeNameNoPV=sim2.nodeName;
%% Max/Min Voltage
f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
h = plot(xrange, Volt_MaxMin_PV(:,1), 'r', xrange, Volt_MaxMin_NoPV(:,1),...
    'g', xrange, Volt_MaxMin_PV(:,2), 'b', xrange, Volt_MaxMin_NoPV(:,2), 'k');
legend(h,{'Max (PV)', 'Max (No PV)','Min (PV)', 'Min (No PV)'},'fontsize',15,'Location','Best');
title('Maximum and Minimum Bus Voltages','fontsize',25)
xlabel('Time of Day','fontsize',25)
ylabel('Voltages (pu)','fontsize',25)
grid on;box on;
xlim([min(xrange) max(xrange)])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

%% Total Loss
% f = figure('units','normalized','outerposition',[0 0 1 1])('visible','off'); hold on;
f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
h = plot(xrange, LossTotal_PV(:,1), 'r', xrange, LossTotal_NoPV(:,1),...
    'g', xrange, LossTotal_PV(:,2), 'b', xrange, LossTotal_NoPV(:,2), 'k');
legend(h,{'kW (PV)', 'kW (No PV)', 'kVar (PV)', 'kVar (No PV)'},'fontsize',15,'Location','Best')
title('Total Losses','fontsize',25)
xlabel('Time of Day','fontsize',25)
xlim([min(xrange) max(xrange)])
ylabel('Losses','fontsize',25)
grid on;box on;
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

%% Line Loss
% f = figure('visible','off'); hold on;
f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
h = plot(xrange, LossLine_PV(:,1), 'r', xrange, LossLine_NoPV(:,1),...
    'g', xrange, LossLine_PV(:,2), 'b', xrange, LossLine_NoPV(:,2), 'k');
legend(h,{'kW (PV)', 'kW (No PV)', 'kVar (PV)', 'kVar (No PV)'},'fontsize',15,'Location','Best')
title('Total Line Losses','fontsize',25)
xlabel('Time of Day','fontsize',25)
xlim([min(xrange) max(xrange)])
ylabel('Losses','fontsize',25)
grid on;box on;
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

%% Total Power
% f = figure('visible','off'); hold on;
if ~isempty(Event_PV.RegControl)
    Names_Reg = fieldnames(Event_PV.RegControl);
    for i=1:length(Names_Reg)
        RegCtr = Names_Reg{i};
        AA.(RegCtr)=unique([Event_PV.RegControl.(RegCtr).Regulations.time]);
    end
end
f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
h = plot(xrange, TotalPower_PV(:,1), 'r', xrange, TotalPower_NoPV(:,1),...
    'g', xrange, TotalPower_PV(:,2), 'b', xrange, TotalPower_NoPV(:,2), 'k');
% % if ~isempty(Event_PV.RegControl)
% %     for ij = 1:length(Names_Reg)
% %         RegCtr = Names_Reg{ij};
% % %         for ii=2:length(AA.(RegCtr))
% % %             line(ones(1,2)*AA.(RegCtr)(ii)/3600,[min(min(min([TotalPower_PV(:,1) ,...
% % %                 TotalPower_PV(:,2)]),min([TotalPower_PV(:,2)]))) max(max([TotalPower_PV(:,1),TotalPower_PV(:,2)]))],...
% % %                 'Color',ColorLine{ij},'LineStyle','--');
% % %             if length(AA.(RegCtr))>2
% % %                 if ii==3
% % %                 text(AA.(RegCtr)(ii)/3600,min(min(min([TotalPower_PV(:,1)]),min([TotalPower_PV(:,2)]))),['\leftarrow Tap op. ' RegCtr],'FontSize',16)
% % %                 end
% % %             end
% % %         end
% %     end
% % end
    legend(h,{'MW (PV)', 'MW (No PV)', 'MVar (PV)', 'MVar (No PV)'},'fontsize',15,'Location','Best')
    title('Total Power','fontsize',25)
    xlabel('Time of Day','fontsize',25)
    xlim([min(xrange) max(xrange)])
    ylabel('Power','fontsize',25)
    grid on;box on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

close all;

f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
h = plot(xrange, -TotalPower_PV(:,1)+TotalPower_NoPV(:,1),...
    'g', xrange, -TotalPower_PV(:,2)+ TotalPower_NoPV(:,2), 'k');
legend(h,{'MW PV prod', 'difference in MVar'},'fontsize',15,'Location','Best')
    title('PV prod','fontsize',25)
    xlabel('Time of Day','fontsize',25)
    xlim([min(xrange) max(xrange)])
    ylabel('Power','fontsize',25)
    grid on;box on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

% % % find transformer at substation in monPV.m1:
% % tt=circuit.transformer(:).sub;
% % for kk=1:length(tt)
% %     ts{kk}=tt{kk}(1);
% % end
% % trsub = find(ismember(ts,'y'));
% % mosub = find(ismember({monPV.m1(:).Name}, ['mon_1_transformer_' circuit.transformer(trsub).Name]));
%% Reg and Cap control 
% Maximum 10 cap or reg curve on one graph but it can be changed.
% f = figure('visible','off'); hold on;
if ~isempty(Event_PV.RegControl)
    numOfRegCurve = min(length( Reg_PV(1,:))-1,length(colorCurve)); 

    romanNum = {'I','II','III','IV','V','VI','VII','VIII','IX','X'};
    LegendList = {};
    for k = 1:numOfRegCurve
        LegendList{k}= {[romanNum{k} '(' num2str(Reg_PV(1,k)) ')']};
    end

    try 
        f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        for z=1:numOfRegCurve
            plot(xrange, Reg_PV((2:end),z),colorCurve{z});
            grid on;box on;
        end
    catch
        for z=1:numOfRegCurve
            plot(xrange, Reg_PV((2:length(Reg_PV)/length(xrange):end),z),colorCurve{z});
            grid on;box on;
        end
    end
    legend([LegendList{:}],'Location','Best')
    title('Voltage Regulator Events (With PV)','fontsize',25)
    xlabel('Time of Day','fontsize',25)
    xlim([min(xrange) max(xrange)])
    ylabel('# Events','fontsize',25)
    grid on;box on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
end

if ~isempty(Event_NoPV.RegControl)
    % f = figure('visible','off'); hold on;
    numOfRegCurve = min(length( Reg_PV(1,:))-1,length(colorCurve)); 

    romanNum = {'I','II','III','IV','V','VI','VII','VIII','IX','X'};
    LegendList = {};
    for k = 1:numOfRegCurve
        LegendList{k}= {[romanNum{k} '(' num2str(Reg_NoPV(1,k)) ')']};
    end

    try 
        f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        for z=1:numOfRegCurve
            plot(xrange, Reg_NoPV((2:end),z),colorCurve{z});
            grid on;box on;
        end
    catch
        for z=1:numOfRegCurve
            plot(xrange, Reg_NoPV((2:length(Reg_NoPV)/length(xrange):end),z),colorCurve{z});
            grid on;box on;
        end
    end
    legend([LegendList{:}],'Location','Best')
    title('Voltage Regulator Events (No PV)','fontsize',25)
    xlabel('Time of Day','fontsize',25)
    xlim([min(xrange) max(xrange)])
    ylabel('# Events','fontsize',25)
    grid on;box on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
end

if ~isempty(Cap_PV)
% f = figure('visible','off'); hold on;
    numOfCapCurve = min(length( Cap_PV(1,:))-1,length(colorCurve)); 

    romanNum = {'I','II','III','IV','V','VI','VII','VIII','IX','X'};
    LegendList = {};
    for k = 1:numOfCapCurve
        LegendList{k}= {[romanNum{k} '(' num2str(Cap_PV(1,k)) ')']};
    end

    f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    for z=1:numOfCapCurve
        plot(xrange, Cap_PV((2:end),z),colorCurve{z});
        grid on;box on;
    end
    legend([LegendList{:}],'Location','Best')
    title('Capacitor Bank Events (With PV)','fontsize',25)
    xlabel('Time of Day','fontsize',25)
    xlim([min(xrange) max(xrange)])
    ylabel('# Events','fontsize',25)
    grid on;box on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
end

if ~isempty(Cap_NoPV)
    % f = figure('visible','off'); hold on;
    numOfCapCurve = min(length( Cap_NoPV(1,:))-1,length(colorCurve)); 

    romanNum = {'I','II','III','IV','V','VI','VII','VIII','IX','X'};
    LegendList = {};
    for k = 1:numOfCapCurve
        LegendList{k}= {[romanNum{k} '(' num2str(Cap_NoPV(1,k)) ')']};
    end

    f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    for z=1:numOfCapCurve
        plot(xrange, Cap_NoPV((2:end),z),colorCurve{z});
        grid on;box on;
    end
    legend([LegendList{:}],'Location','Best')
    title('Capacitor Bank Events (No PV)','fontsize',25)
    xlabel('Time of Day','fontsize',25)
    xlim([min(xrange) max(xrange)])
    ylabel('# Events','fontsize',25)
    grid on;box on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
end

if ~isempty(Reg_PV)
    RegTotal_PV = sum(Reg_PV(2:end,end));
    RegTotal_NoPV = sum(Reg_NoPV(2:end,end));
    % f = figure('visible','off'); hold on;
    try 
        f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        h = plot(xrange, Reg_PV(2:end,end), 'r', xrange, Reg_NoPV(2:end,end), 'g');
    catch
        h = plot(xrange, Reg_PV(2:length(Reg_PV)/length(xrange):length(Reg_PV),end), 'r', xrange, Reg_NoPV(2:length(Reg_NoPV)/length(xrange):length(Reg_NoPV),end), 'g');
    end
    legend(h,{['Total_P_V (' num2str(RegTotal_PV) ')'], ['Total_N_o_P_V (' num2str(RegTotal_NoPV) ')']},'Location','Best')
    title('Voltage Regulator Events - Total','fontsize',25)
    xlabel('Time of Day','fontsize',25)
    xlim([min(xrange) max(xrange)])
    ylabel('# Events','fontsize',25)
    grid on;box on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
end
if ~isempty(Cap_PV)
    CapTotal_PV = sum(Cap_PV(2:end,end));
    CapTotal_NoPV = sum(Cap_NoPV(2:end,end));
    % f = figure('visible','off'); hold on;
    f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    h = plot(xrange, Cap_PV(2:end,end), 'r', xrange, Cap_NoPV(2:end,end), 'g');
    legend(h,{['Total_P_V (' num2str(CapTotal_PV) ')'], ['Total_N_o_P_V (' num2str(CapTotal_NoPV) ')']},'Location','Best')
    title('Capacitor Bank Events - Total','fontsize',25)
    xlabel('Time of Day','fontsize',25)
    xlim([min(xrange) max(xrange)])
    ylabel('# Events','fontsize',25)
    grid on;box on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
end
close all;
%% Tap operation plots
if ~isempty(tap2plotPV)
    RegName = fieldnames(tap2plotPV);
    for nbReg=1:length(RegName)
        f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        h = plot(xrange,[tap2plotPV.(RegName{nbReg})],'b',xrange,[tap2plotNoPV.(RegName{nbReg})],'r');
        legend(h,{['With PV'],['Without PV']},'Location','Best')
        title(['Voltage Regulator ' num2str(nbReg) ' (' RegName{nbReg} ') Tap positions'],'fontsize',25)
        xlabel('Time of Day','fontsize',25)
        xlim([min(xrange) max(xrange)])
        ylabel('(pu)','fontsize',25)
        grid on;box on;
        saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
        saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
    end
    
    romanNum = {'I','II','III','IV','V','VI','VII','VIII','IX','X'};
    LegendList = {};
    for k = 1:min(length(RegName),10)
        LegendList{k}= {[romanNum{k} ' : ' RegName{k}]};
    end

    f = figure('units','normalized','outerposition',[0 0 1 1]); 
    subplot(2,1,2,'align'), 
    hold on;
    for z=1:length(RegName)
        plot(xrange,[tap2plotPV.(RegName{z})],colorCurve{z});
        grid on;box on;
    end
    legend([LegendList{:}],'fontsize',25)
    grid on;box on;
    xlabel('Time of Day','fontsize',25)
    xlim([min(xrange) max(xrange)])
    ylabel('Voltage Regulators Tap positions (pu)','fontsize',25)
    subplot(2,1,1,'align'), plot(xrange, TotalPower_PV(:,1), 'r', xrange, TotalPower_NoPV(:,1),...
        'g', xrange, TotalPower_PV(:,2), 'b', xrange, TotalPower_NoPV(:,2), 'k');
    legend({'MW (PV)', 'MW (No PV)', 'MVar (PV)', 'MVar (No PV)'},'fontsize',15)
    title('Total Power compare to tap position','fontsize',25)
    xlim([min(xrange) max(xrange)])
    grid on;box on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
end
close all;
%% Cap step plot
if ~isempty(cap2plotPV)
    CapName = fieldnames(cap2plotPV);
    for nbCap=1:length(CapName)
        f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        h = plot(xrange,[cap2plotPV.(CapName{nbCap})],'b',xrange,[cap2plotNoPV.(CapName{nbCap})],'r');
        legend(h,{['With PV'],['Without PV']},'Location','Best')
        title(['Capacitor Bank ' num2str(nbCap) ' (' CapName{nbCap} ') step Positions'],'fontsize',25)
        xlabel('Time of Day','fontsize',25)
        xlim([min(xrange) max(xrange)])
        ylabel('Step','fontsize',25)
        grid on;box on;
        saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
        saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
    end
end

close all;
%% Plot voltage along the circuit
% change it to plot or not voltage along the feeder at each hour. *****
if exist('plotVoltageProfile','var') && plotVoltageProfile
    volt_pv_plot = Volt_PV;
    volt_pv_plot([2<Volt_PV]) = nan;
    volt_pv_plot([0.5>Volt_PV]) = nan;
    volt_nopv_plot = Volt_NoPV;
    volt_nopv_plot([2<Volt_NoPV]) = nan;
    volt_nopv_plot([0.5>Volt_NoPV]) = nan;
    %% pb voltage is ..x.. 
%     if length(DistPV)>length(DistNoPV)
    if length(volt_pv_plot)<length(volt_nopv_plot)
%         dist_plot = DistNoPV(~isnan(volt_pv_plot(1,:))); % DistPV and DistNoPV 
%         may not have the same length. Don't know why! So to make it work
%         I use the shortest one. It adds errors but otherwise it doesn't
%         work!
		[aa bb]=ismember(NodeNamePV,NodeNameNoPV);
		volt_nopv_plot=volt_nopv_plot(:,bb);
		dist_plot = DistPV(aa);
    else
        [aa bb]=ismember(NodeNameNoPV,NodeNamePV);
        volt_pv_plot=volt_pv_plot(:,bb);
        dist_plot = DistNoPV(aa);
    end
%     end
    
    [order order] = sort(dist_plot);
    
    %% Plot voltage 
    if ~exist('plotTime','var')
        timearray= 120:120:2880;
    else 
        timearray=plotTime*120;
    end
    for time= timearray;
        f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        for i=1:length(dist_plot)
            plot([dist_plot(order(i)) dist_plot(order(i))],[volt_pv_plot(time,order(i)) volt_nopv_plot(time,order(i))] ,'r');
        end
        h = plot(dist_plot(order),volt_pv_plot(time,order),'x',dist_plot(order),volt_nopv_plot(time,order),'.');
        legend(h,{'With PV','Without PV'},'fontsize',15);
        grid on;box on;
        xlim([0 max(dist_plot)]);ylim([min(min(volt_pv_plot(volt_pv_plot>0.4)),min(volt_nopv_plot(volt_nopv_plot>0.4)))  max(max(max(volt_pv_plot)),max(max(volt_nopv_plot)))])
        xlabel('Distance, km','fontsize',20);
        ylabel('Voltage, pu','fontsize',20);
        title(['Voltage at ' num2str(time/120)],'fontsize',20)
        saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
        saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
        if (time==120*8) || (time==120*16)
             close all;
        end
    end    
end
