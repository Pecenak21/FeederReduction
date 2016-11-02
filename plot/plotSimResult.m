function [result_dir,d] = plotSimResult( res, type, plotTime, result_dir, overwrite, ordering, substation )
%  plotSimResult Function to plot results of dssSimulation
% INPUT:
%   - res: result struct obtained from running simMain
%   - result_dir: Path to where to save figures
%   - type:   * single = Plots for a single simulation result struct
%                        The time series plots include min/max Volt, lineloss, total loss,
%                        power consumption, convergence, and transformers'
%                        and capacitors' taps position during the time of simulation.
%             * multiple = plot for a result struct (from simMain) with multiple simulations
%                        The plots are sum plots (add all values for the whole time of simulation)
%                        of all values in 'single' plots versus different PV penetrations
%                        (i.e. x-axis is PV pen level) for each feeder.
%             * voltageProfile = plot voltage profile for all simulation results for time
%                        specified in plotTime parameter. plotTime can be an array of time.
%   - plotTime: (optionnal) time at which you want to plot.
%               can be an array with different time or time range.

%% initialize
conf = getConf;
if ~exist('result_dir','var') || isempty(result_dir)
    result_dir = [conf.outputDir '/fig'];
    if ~exist(result_dir,'dir'), mkdir(result_dir); end
end
if ~exist('overwrite','var') || isempty(overwrite), overwrite = 1; end
if ~exist('type','var') || isempty(type)
    type = {'single','multiple'};
elseif exist('type','var') && ~isempty(type) && ~iscell(type)
    type = {type};
end
type = lower(type);

t=res.time; % time
fnPrefix = fNamePrefix('',res,result_dir);

%% intialize output
d = [];

%% find 0% penetration index if exist
if length(res) > 1 && ismember(0,[res.penLevel])
    [~,idPen0] = ismember(0,[res.penLevel]);
end

%% plots for each simulation results
if ismember('single',type)
    f = figure(10);
    for i = 1:length(res)
        r = res(i);
        %% Max/Min Voltage
        fn = [fnPrefix{i} '_VoltMinMax'];
        if overwrite || ~exist([fn '.png'],'file') || ~exist([fn '.fig'],'file')
            figure(10); clf; h=plot(t, r.VoltMaxMin(:,1), t, r.VoltMaxMin(:,2)); datetick; grid on;box on;
            legend(h,{'Max', 'Min'},'Location','Best');xlabel('Time, [HH:MM]'); ylabel('Voltage, [p.u.]')
            
            % save
            saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
        end
        
        %% Total Loss
        fn = [fnPrefix{i} '_TotalLoss'];
        if overwrite || ~exist([fn '.png'],'file') || ~exist([fn '.fig'],'file')
            figure(10); clf; h = plot(t, r.TotalLoss(:,1)/1000, t, r.TotalLoss(:,2)/1000);datetick; grid on;box on;
            legend(h,{'kW','kVar'},'Location','Best'); xlabel('Time, [HH:MM]'); ylabel('Total Losses, [kX]')
            
            % save
            saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
        end
        
        %% Line Losses
        fn = [fnPrefix{i} '_LineLoss'];
        if overwrite || ~exist([fn '.png'],'file') || ~exist([fn '.fig'],'file')
            figure(10); clf; h = plot(t, r.LineLoss(:,1)/1000, t, r.LineLoss(:,2)/1000); datetick; grid on;box on;
            legend(h,{'kW','kVar'},'Location','Best'); xlabel('Time, [HH:MM]'), ylabel('Line Losses, [kX]')
            
            % save
            saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
        end
        
        %% convergence
        fn = [fnPrefix{i} '_Convergence'];
        if overwrite || ~exist([fn '.png'],'file') || ~exist([fn '.fig'],'file')
            figure(10); clf; h = plot(t, r.converged,'-*'); datetick; grid on;box on;
            legend(h,{'0-Not converged; 1-Converged; 2-Not controlly converged'},'Location','Best')
            xlabel('Time, [HH:MM]'), ylabel('Convergence Index, [-]')
            
            saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
        end
        
        %% Total Power
        fn = [fnPrefix{i} '_NetPower'];
        if overwrite || ~exist([fn '.png'],'file') || ~exist([fn '.fig'],'file')
            figure(10); clf; h = plot(t, r.TotalPower(:,1), t, r.TotalPower(:,2)); datetick; grid on;box on;
            legend(h,{'MW', 'MVar'},'Location','Best')
            xlabel('Time, [HH:MM]'), ylabel('Net Power, [MX]')
            
            saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
        end
        
        %% Reg and Cap control
        if isfield(res,'tapPos')
            f1=figure(101); f2=figure(102);
            for j = 1:length(r.tapPos.transformer)
                fn1 = sprintf('%s_%s%03d_%s',fnPrefix{i},'Tranx',j,'TapPos');
                fn2 = sprintf('%s_%s%03d_%s',fnPrefix{i},'Tranx',j,'Volt');
                if overwrite || ~exist([fn1 '.png'],'file') || ~exist([fn1 '.fig'],'file') || ~exist([fn2 '.png'],'file') || ~exist([fn2 '.fig'],'file')
                    % tranx tap position
                    figure(101); clf; h1 = plot(r.time,r.tapPos.transformer(j).pos);datetick; grid on;box on; hold off;
                    legend(h1,{'Winding 1','Winding 2','Winding 3'},'Location','Best');
                    xlabel('Time, [HH:MM]'), ylabel(['Transformer ' num2str(j) '''s Tap Position, [-]']);
                    % save
                    
                    saveFigure(f1,[fn1 '.png'],overwrite); saveFigure(f1,[fn1 '.fig'],overwrite);
                    
                    % tranx voltage
                    figure(102); clf; h2 = plot(r.time,r.tapPos.transformer(j).V);datetick; grid on;box on; hold off;
                    legend(h2,{'Winding 1','Winding 2','Winding 3'},'Location','Best')
                    xlabel('Time, [HH:MM]'), ylabel(['Transformer ' num2str(j) '''s Voltage, [p.u.]']);
                    % save
                    saveFigure(f2,[fn2 '.png'],overwrite); saveFigure(f2,[fn2 '.fig'],overwrite);
                end
            end
        end
        
        %% capacitor taps
        if isfield(r,'tapPos') && isfield(r.tapPos,'capacitor')
        end
    end
end

%% PV production for all except the 0% pv penetration result
if ismember('pvproduction',type) && length(res)>1 && exist('idPen0','var')
%     id = setdiff(1:length(res),idPen0);
%     for i=1:length(id)
%         idx = id(i);
%         f = figure(10); grid on; box on;
%         h = plot(t, -res(idx).TotalPower(:,1)+res(idPen0).TotalPower(:,1),t, -res(idx).TotalPower(:,2)+res(idPen0).TotalPower(:,2)); datetick;
%         legend(h,{'MW', 'MVar'},'Location','Best')
%         xlabel('Time, [HH:MM]'), ylabel('PV Production, [MX]')
%
%         fn = [fnPrefix{i} '_PVProduction'];saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
%     end
end

%% plot group plots
if ismember('multiple',type) && length(res) > 1
    if exist('plotTime','var') && ~isempty(plotTime)
        tid = (res(1).time >= min(plotTime)) & (res(1).time <= max(plotTime));
        % if tid is empty then there might be mismatch in day, handle this by removing the day dependence and focus on hours of interest
        if isempty(find(tid,1))
            plotTime = plotTime - floor(plotTime(1));
            t1 = res(1).time; t1 = t1 - floor(t1(1));
            tid = (t1 >= min(plotTime)) & (t1 <= max(plotTime));
        end
    end
    % process data
    % x-axis: pv penetration
    pen = sort(unique([res.penLevel]));
    % create unique id
    % count number of unique (feeder + feederSetup)
    ufd = cell(length(res),1);
    for i = 1:length(res)
        ufd{i} = lower([res(i).feederName '_' res(i).feederSetup '_' res(i).timeDay]);
    end
    ufd = unique(ufd);
    d.totLoadKw = nan(length(ufd),length(pen));
    d.totLoadKvar = nan(length(ufd),length(pen));
    d.totTapOpe = nan(length(ufd),length(pen));
    d.totLineLossVA = nan(length(ufd),length(pen));
    d.totLossVA = nan(length(ufd),length(pen));
    d.Vmin = nan(length(ufd),length(pen));
    d.VminTime = nan(length(ufd),length(pen));
    d.Vmax = nan(length(ufd),length(pen));
    d.VmaxTime = nan(length(ufd),length(pen));
    d.VsubPhase1Std = nan(length(ufd),length(pen));
    d.VsubPhase1Std = nan(length(ufd),length(pen)); 
    d.VsubPhase1DiffStd = nan(length(ufd),length(pen));
    
    % create index for the 0% penetration
    id0 = nan(length(ufd),1);
    for i = 1:length(res)
        r = res(i);
        [~,id1] = ismember(lower([res(i).feederName '_' res(i).feederSetup '_' res(i).timeDay]),ufd);
        if isequal(r.penLevel,0)
            id0(id1) = i;
        end
    end
    for i = 1:length(res)
        r = res(i);
        [~,id1] = ismember(lower([res(i).feederName '_' res(i).feederSetup '_' res(i).timeDay]),ufd);
        [~,id2] = ismember(r.penLevel,pen);
        if exist('tid','var')
            t = res(i).time(tid);
            d.totLoadKw(id1,id2) = sum(r.TotalPower(tid,1));
            d.totLoadKvar(id1,id2) = sum(r.TotalPower(tid,2));
            d.totTapOpe(id1,id2) = r.totTranxTapOpe + r.totCapTapOpe;
            [d.Vmax(id1,id2), id] = max(r.VoltMaxMin(tid,1));
            d.VmaxTime(id1,id2) = t(id);
            [d.Vmin(id1,id2), id] = min(r.VoltMaxMin(tid,2));
            d.VminTime(id1,id2) = t(id);
            d.totLossVA(id1,id2) = sum(r.TotalLoss(tid,1));
            d.totLineLossVA(id1,id2) = sum(r.LineLoss(tid,1));
            d.VsubPhase1Std(id1,id2) = std(r.Voltage(tid,4)); % use the voltage of the 4th node 
            d.VsubPhase1DiffStd(id1,id2) = std(abs(r.Voltage(tid,4)-res(id0(id1)).Voltage(tid,4)));
            
        else 
            d.totLoadKw(id1,id2) = sum(r.TotalPower(:,1));
            d.totLoadKvar(id1,id2) = sum(r.TotalPower(:,2));
            d.totTapOpe(id1,id2) = r.totTranxTapOpe + r.totCapTapOpe;
            [d.Vmax(id1,id2), id] = max(r.VoltMaxMin(:,1));
            d.VmaxTime(id1,id2) = res(i).time(id);
            [d.Vmin(id1,id2), id] = min(r.VoltMaxMin(:,2));
            d.VminTime(id1,id2) = res(i).time(id);
            d.totLossVA(id1,id2) = sum(r.TotalLoss(:,1));
            d.totLineLossVA(id1,id2) = sum(r.LineLoss(:,1));
            d.VsubPhase1Std(id1,id2) = std(r.Voltage(:,4));
            d.VsubPhase1DiffStd(id1,id2) = std(abs( r.Voltage(:,4 )- res(id0(id1)).Voltage(:,4) ));
        end
    end
    if exist('ordering','var') && ~isempty(ordering)
        fns = fieldnames(d);
        for i = 1:length(fns)
            fn = fns{i};
            d.(fn) = d.(fn)(ordering,:);
        end
        ufd = ufd(ordering);
    end
    %% daily Load consumption (kw)
    f21=figure(21); %set(gca,'fontsize',12);
    h = bar(pen, d.totLoadKw','group','basevalue',0); hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    legend(ufd,'location','Best','interpreter','none'); box on;set(gca,'Ygrid','on');
    %xlim([-15 115]);
    xlabel('PV Penetration, [%]'); ylabel('Total Net Consumption, [MWh]');
    set(gcf,'color','w')
    set(gcf,'position', [420   270   440   309]); %ylim([0 510])
    fn = [result_dir '/PVpen_totNetPowConsMW'];
    if ~exist(result_dir,'dir'), mkdir(result_dir); end
    saveFigure(f21,[fn '.png'],overwrite); saveFigure(f21,[fn '.fig'],overwrite);
    
    %% total consumption saving in percent and bar graph style
    f22=figure(22); %set(gca,'fontsize',12)
    x = d.totLoadKw ; x = 100 - x./ repmat(x(:,1),[1 size(x,2)])*100;
    h = bar(pen, x','group','basevalue',0); hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    legend(ufd,'location','Best','interpreter','none'); box on;set(gca,'Ygrid','on');
    %xlim([-15 115])
    %set(gca,'XTick',0:25:100);
    xlabel('PV Penetration, [%]'); ylabel('Consumption Saving, [%]'); set(gcf,'color','w')
    set(gcf,'position', [420   270   440   309]);% ylim([0 120])
    fn = [result_dir '/PVpen_totNetPowConsPercentBar'];
    saveFigure(f22,[fn '.png'],overwrite); saveFigure(f22,[fn '.fig'],overwrite);
    
    %% daily consumption saving in percent, line chart
    f23=figure(23);
    x = d.totLoadKw ; x = 100 - x./ repmat(x(:,1),[1 size(x,2)])*100;
    plot(pen,x','*-','linewidth',2);
    legend(ufd,'location','Best','interpreter','none'); box on; grid on;
    xlabel('PV Penetration, [%]'); ylabel('Consumption Saving, [%]');
    set(gcf,'position', [420   270   440   309]); set(gcf,'color','w'); %ylim([0 19])
    %set(gca,'XTick',0:25:100);
    fn = [result_dir '/PVpen_totNetPowConsPercentLine'];
    saveFigure(f23,[fn '.png'],overwrite); saveFigure(f23,[fn '.fig'],overwrite);
    
    %% Voltage max/min figure
    % bar chart
    f24 = figure(24); hold off;
    h = bar(pen, d.Vmax','group','basevalue',1,'barwidth',.7); hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    % add time on top of the bars
    for i1=1:length(h)
        for i2 = 1:length(h(i1).XData)
            text(h(i1).XData(i2) + h(i1).XOffset, h(i1).YData(i2), datestr(d.VmaxTime(i1,i2),' HH:MM'),...
                'HorizontalAlignment','left',...
                'VerticalAlignment','middle','rotation',90,'fontsize',8)
        end
    end
    hold on;
    h = bar(pen,d.Vmin','group','basevalue',1,'barwidth',.7);
    hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    legend(ufd,'location','Best','interpreter','none'); box on; grid on;
    % add time to the bottom of the bars
    for i1=1:length(h)
        for i2 = 1:length(h(i1).XData)
            text(h(i1).XData(i2) + h(i1).XOffset, h(i1).YData(i2), datestr(d.VminTime(i1,i2),'HH:MM '),...
                'HorizontalAlignment','right',...
                'VerticalAlignment','middle','rotation',90,'fontsize',8)
        end
    end
    %xlim([-15 115])
    %set(gca,'XTick',0:25:100);
    xlabel('PV Penetration, [%]'); ylabel('Voltage, [p.u.]'); set(gcf,'color','w')
    set(gcf,'position', [420   270   440   309]);
    yl = ylim; ylim([yl(1)-0.02, yl(2)+.02]);
    fn = [result_dir '/PVpen_VMaxMin'];
    saveFigure(f24,[fn '.png'],overwrite); saveFigure(f24,[fn '.fig'],overwrite);
    hold off;
    
    %% voltage max - line
    f = figure; 
    plot(pen,d.Vmax','*-','linewidth',1.5);
    resetFigure(f);
    legend(ufd,'location','Best','interpreter','none'); box on; grid on;
    xlabel('PV penetration, [%]'); ylabel('Max Voltage, [p.u.]');
    fn = [result_dir '/PVpen_VMax_Line'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    hold off;
    
    %% voltage max - bar
    f = figure; 
    bar(pen,d.Vmax','group','basevalue',1);
    resetFigure(f);
    legend(ufd,'location','Best','interpreter','none'); box on; grid on;
    xlabel('PV penetration, [%]'); ylabel('Max Voltage, [p.u.]');
    fn = [result_dir '/PVpen_VMax_Bar'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    hold off;
    
    %% voltage min - line
    f = figure; 
    plot(pen,d.Vmin','*-','linewidth',1.5);
    resetFigure(f);
    legend(ufd,'location','Best','interpreter','none'); box on; grid on;
    xlabel('PV penetration, [%]'); ylabel('Min Voltage, [p.u.]');
    fn = [result_dir '/PVpen_VMin_Line'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    hold off;
    
    %% voltage min - bar
    f = figure; 
    bar(pen,d.Vmin','group','basevalue',1);
    resetFigure(f);
    legend(ufd,'location','Best','interpreter','none'); box on; grid on;
    xlabel('PV penetration, [%]'); ylabel('Min Voltage, [p.u.]');
    fn = [result_dir '/PVpen_VMin_Bar'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    hold off;
    
    %% tap operations: regulator and capacitor operations
    f25 = figure(25); 
    h = bar(pen, d.totTapOpe','group','basevalue',0); hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    legend(ufd,'location','Best','interpreter','none'); box on;set(gca,'Ygrid','on');
    %xlim([-15 115]); set(gca,'XTick',0:25:100);
    xlabel('PV Penetration, [%]'); ylabel('Tap Operations, [-]'); set(gcf,'color','w')
    set(gcf,'position', [420   270   440   309]); %ylim([0 1050])
    fn = [result_dir '/PVpen_tapOperNum'];
    saveFigure(f25,[fn '.png'],overwrite); saveFigure(f25,[fn '.fig'],overwrite);
    
    % tap operation [%]
    f26 = figure(26); set(gca,'fontsize',12);
    x = d.totTapOpe; x = x./ repmat(x(:,1),[1 size(x,2)])*100;
    h = bar(pen, x','group','basevalue',0); hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    legend(ufd,'location','Best','interpreter','none'); box on;set(gca,'Ygrid','on');
    %xlim([-15 115]); set(gca,'XTick',0:25:100);
    xlabel('PV Penetration, [%]'); ylabel('Tap Operations, [%]'); set(gcf,'color','w')
    set(gcf,'position', [420   270   440   309]);% ylim([0 400])
    fn = [result_dir '/PVpen_tapOperPercentBar'];
    saveFigure(f26,[fn '.png'],overwrite); saveFigure(f26,[fn '.fig'],overwrite);
    
    %% tap operation [%], line chart
    f27 = figure(27);
    x = d.totTapOpe;
    % for cases with 0 tap operations, add 1 to all PV penetration 
    oId = sum((d.totTapOpe==0),2)>0;
    if sum(oId)>0
        x(oId,:) = x(oId,:) + 1;
    end
    x = x./ repmat(x(:,1),[1 size(x,2)])*100;
    plot(pen,x','*-','linewidth',2);
    legend(ufd,'location','Best','interpreter','none'); grid on; box on;
    xlabel('PV Penetration, [%]'); ylabel('Tap Operations, [%]');
    set(gcf,'position', [420   270   440   309]); set(gcf,'color','w'); %ylim([30 400])
    %set(gca,'XTick',0:25:100);
    fn = [result_dir '/PVpen_tapOperPercentLine'];
    saveFigure(f27,[fn '.png'],overwrite); saveFigure(f27,[fn '.fig'],overwrite);
    
    %% tap operation normalized [-], line chart
    f = figure(127);
    x = d.totTapOpe;
    % for cases with 0 tap operations, add 1 to all PV penetration 
    oId = sum((d.totTapOpe==0),2)>0;
    if sum(oId)>0
        x(oId,:) = x(oId,:) + 1;
    end
    x = x./ repmat(x(:,1),[1 size(x,2)]);
    plot(pen,x','*-','linewidth',2);
    legend(ufd,'location','Best','interpreter','none'); grid on; box on;
    xlabel('PV Penetration, [%]'); ylabel('Normalized Tap Operations, [-]');
    set(gcf,'position', [420   270   440   309]); set(gcf,'color','w'); %ylim([30 400])
    %set(gca,'XTick',0:25:100);
    fn = [result_dir '/PVpen_tapOperNormalizedLine'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    
    %% tap operation normalized [-], bar chart
    f = figure;
    x = d.totTapOpe;
    % for cases with 0 tap operations, add 1 to all PV penetration 
    oId = sum((d.totTapOpe==0),2)>0;
    if sum(oId)>0
        x(oId,:) = x(oId,:) + 1;
    end
    x = x./ repmat(x(:,1),[1 size(x,2)]);
    bar(pen,x','group','basevalue',0); hBaseline = get(h,'baseline');
    resetFigure(f); 
    legend(ufd,'location','Best','interpreter','none'); grid on; box on;
    xlabel('PV Penetration, [%]'); ylabel('Normalized Tap Operations, [-]');
    fn = [result_dir '/PVpen_tapOperNormalizedBar'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    
    %% total loss in bar style
    f28 = figure(28); h = bar(pen, d.totLossVA'/10^6,'group','basevalue',0); hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    legend(ufd,'location','Best','interpreter','none'); box on;set(gca,'Ygrid','on');
    %xlim([-15 115]); set(gca,'XTick',0:25:100);
    xlabel('PV Penetration, [%]'); ylabel('Total Losses, [MWh]'); set(gcf,'color','w')
    set(gcf,'position', [420   270   440   309]); %ylim([0 9])
    fn = [result_dir '/PVpen_totLossBar'];
    saveFigure(f28,[fn '.png'],overwrite); saveFigure(f28,[fn '.fig'],overwrite);
    
    %% total loss reduction in percent
    f29 = figure(29); set(gca,'fontsize',12)
    x = d.totLossVA/10^6; x = 100 - x./ repmat(x(:,1),[1 size(x,2)])*100;
    h = bar(pen, x','group','basevalue',0); hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    legend(ufd,'location','Best','interpreter','none'); box on;set(gca,'Ygrid','on');
    
    %xlim([-15 115])
    set(gca,'XTick',0:25:100);
    xlabel('PV Penetration, [%]'); ylabel('Total Loss Reduction, [%]'); set(gcf,'color','w')
    set(gcf,'position', [420   270   440   309]); %ylim([0 17])
    fn = [result_dir '/PVpen_totLossReduc'];
    saveFigure(f29,[fn '.png'],overwrite); saveFigure(f29,[fn '.fig'],overwrite);
    
    %% total loss reduction in percent, line chart
    x = d.totLossVA/10^6; x = 100 - x./ repmat(x(:,1),[1 size(x,2)])*100;
    f30 = figure(30); plot(pen,x','*-','linewidth',2);
    legend(ufd,'location','Best','interpreter','none'); grid on; box on;
    xlabel('PV Penetration, [%]'); ylabel('Total Loss Reduction, [%]');
    set(gcf,'position', [420   270   440   309]); set(gcf,'color','w'); %ylim([-5 24])
    set(gca,'XTick',0:25:100);
    fn = [result_dir '/PVpen_totLossReducPercentLine'];
    saveFigure(f30,[fn '.png'],overwrite); saveFigure(f30,[fn '.fig'],overwrite);
    
    %% total loss normalized in percent, line chart
    x = d.totLossVA/10^6; x = x./ repmat(x(:,1),[1 size(x,2)])*100;
    f = figure; plot(pen,x','*-','linewidth',2);
    legend(ufd,'location','Best','interpreter','none'); grid on; box on;
    xlabel('PV Penetration, [%]'); ylabel('Nomarlized Total Losses, [%]');
    set(gcf,'position', [420   270   440   309]); set(gcf,'color','w'); %ylim([-5 24])
    set(gca,'XTick',0:25:100);
    fn = [result_dir '/PVpen_totLossNormalizedPercentLine'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    
    %% lineloss in bar chart
    f31 = figure(31); h = bar(pen, d.totLineLossVA'/10^6,'group','basevalue',0); hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    legend(ufd,'location','Best','interpreter','none'); box on;set(gca,'Ygrid','on');
    %xlim([-15 115]);set(gca,'XTick',0:25:100);
    xlabel('PV Penetration, [%]'); ylabel('Line Losses, [MWh]'); set(gcf,'color','w')
    set(gcf,'position', [420   270   440   309]); %ylim([0 9])
    fn = [result_dir '/PVpen_lineLoss'];
    saveFigure(f31,[fn '.png'],overwrite); saveFigure(f31,[fn '.fig'],overwrite);
    
    %% lineloss reduction in percent
    x = d.totLineLossVA/10^6; x = 100 - x./ repmat(x(:,1),[1 size(x,2)])*100;
    f32 = figure(32); h = bar(pen, x','group','basevalue',0); hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    legend(ufd,'location','Best','interpreter','none'); box on;set(gca,'Ygrid','on');
    %xlim([-15 115]), set(gca,'XTick',0:25:100);
    xlabel('PV Penetration, [%]'); ylabel('Line Loss Reduction, [%]'); set(gcf,'color','w')
    set(gcf,'position', [420   270   440   309]); %ylim([0 21])
    fn = [result_dir '/PVpen_linelossReducPercentBar'];
    saveFigure(f32,[fn '.png'],overwrite); saveFigure(f32,[fn '.fig'],overwrite);
    
    %% lineloss reduction in percent, line chart
    x = d.totLineLossVA/10^6; x = 100 - x./ repmat(x(:,1),[1 size(x,2)])*100;
    f33 = figure(33); plot(pen,x','*-','linewidth',2);
    legend(ufd,'location','Best','interpreter','none'); grid on; box on;
    xlabel('PV Penetration, [%]'); ylabel('Line Loss Reduction, [%]');
    set(gcf,'position', [420   270   440   309]); set(gcf,'color','w'); %ylim([-5 23])
    %set(gca,'XTick',0:25:100);
    fn = [result_dir '/PVpen_linelossReducPercentLine'];
    saveFigure(f33,[fn '.png'],overwrite); saveFigure(f33,[fn '.fig'],overwrite);
    
    %% lineloss normalized in percent, line chart
    x = d.totLineLossVA/10^6; x = x./ repmat(x(:,1),[1 size(x,2)])*100;
    f = figure; plot(pen,x','*-','linewidth',2);
    legend(ufd,'location','Best','interpreter','none'); grid on; box on;
    xlabel('PV Penetration, [%]'); ylabel('Normalized Linelosses, [%]');
    set(gcf,'position', [420   270   440   309]); set(gcf,'color','w'); %ylim([-5 23])
    %set(gca,'XTick',0:25:100);
    fn = [result_dir '/PVpen_linelossNormalizedPercentLine'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    
    %% lineloss normalized in percent, bar chart
    x = d.totLineLossVA/10^6; x = x./ repmat(x(:,1),[1 size(x,2)])*100;
    f = figure; bar(pen,x');
    legend(ufd,'location','Best','interpreter','none'); grid on; box on;
    xlabel('PV Penetration, [%]'); ylabel('Normalized Linelosses, [%]');
    set(gcf,'position', [420   270   440   309]); set(gcf,'color','w'); %ylim([-5 23])
    %set(gca,'XTick',0:25:100);
    fn = [result_dir '/PVpen_linelossNormalizedPercentBar'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    
    %% substation voltage volativity index (standard deviation)
    x = d.VsubPhase1Std;
    f = figure(34); plot(pen,x','*-','linewidth',2);
    legend(ufd,'location','Best','interpreter','none'); grid on; box on;
    xlabel('PV Penetration, [%]'); ylabel('Voltage Volatility, [-]');
    resetFigure(gcf);
    fn = [result_dir '/PVpen_voltageVolatilityStd'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    
    %%
    f = figure(36); h = bar(pen, d.VsubPhase1Std','group','basevalue',0); hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    legend(ufd,'location','Best','interpreter','none'); box on;set(gca,'Ygrid','on');
    %xlim([-15 115]);set(gca,'XTick',0:25:100);
    xlabel('PV Penetration, [%]'); ylabel('Voltage Volatility, [-]');
    resetFigure(gcf);
    fn = [result_dir '/PVpen_voltageVolatilityStdBar'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    
    %% substation voltage volativity index 2 (standard deviation of high PV pen case - 0% PV based case)
    x = d.VsubPhase1DiffStd;
    f = figure(35); plot(pen,x','*-','linewidth',2);
    legend(ufd,'location','Best','interpreter','none'); grid on; box on;
    xlabel('PV Penetration, [%]'); ylabel('Voltage Volatility, [-]');
    set(gcf,'position', [420   270   440   309]); set(gcf,'color','w'); %ylim([-5 23])
    %set(gca,'XTick',0:25:100);
    fn = [result_dir '/PVpen_voltageVolatilityDiffStd'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
    
    %%
    f = figure(37); h = bar(pen, d.VsubPhase1DiffStd','group','basevalue',0); hBaseline = get(h,'baseline');
    for i = 1:length(hBaseline)
        set(hBaseline{i},'visible','off');
    end
    legend(ufd,'location','Best','interpreter','none'); box on;set(gca,'Ygrid','on');
    %xlim([-15 115]);set(gca,'XTick',0:25:100);
    xlabel('PV Penetration, [%]'); ylabel('Voltage Volatility, [-]');
    resetFigure(gcf);
    fn = [result_dir '/PVpen_voltageVolatilityDiffStdBar'];
    saveFigure(f,[fn '.png'],overwrite); saveFigure(f,[fn '.fig'],overwrite);
end
end