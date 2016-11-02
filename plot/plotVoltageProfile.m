function [o, dat, fh] = plotVoltageProfile(o,ignoreBus,h,doplot)
% plot voltage profile with busName as datatip
if ~isa(o,'COM.OpendssEngine_dss')
    o = dssget(o,[],[],[],1);
end

if ~exist('doplot','var')
    doplot = 1;
end
%% Plot voltage
if doplot
    if exist('h','var') && ~isempty(h)
        fh = h;
    else
        fh = figure; 
    end
end
d = cell(1,3); v = cell(1,3); n = cell(1,3);
for i = 1:3
    d{i} = o.ActiveCircuit.AllNodeDistancesByPhase(i);
    v{i} = o.ActiveCircuit.AllNodeVmagPUByPhase(i);
    n{i} = o.ActiveCircuit.AllNodeNamesByPhase(i);
    if exist('ignoreBus','var')
        ignoreBusId = ismember(lower(cleanBus(n{i})),lower(cleanBus(ignoreBus)));
        n{i} = n{i}(~ignoreBusId);
        d{i} = d{i}(~ignoreBusId);
        v{i} = v{i}(~ignoreBusId);
    end
end
dat.nodeName = [n{1};n{2};n{3}];
dat.nodeV = [v{1},v{2},v{3}];
dat.nodeDist = [d{1},d{2},d{3}];
dat.numNodeInEachPhase = [length([v{1}]) length([v{2}]) length([v{3}])];

if doplot
    plot(d{1},v{1},'.',d{2},v{2},'.',d{3},v{3},'.','linewidth',10);
    legend({'Phase 1','Phase 2','Phase 3'});
    grid on;box on;
    % % xlim([0 max(dist_plot)]);ylim([min(min(volt_pv_plot(volt_pv_plot>0.4)),min(volt_nopv_plot(volt_nopv_plot>0.4)))  max(max(max(volt_pv_plot)),max(max(volt_nopv_plot)))])
    xlabel('Distance [km]');
    ylabel('Voltage [p.u.]');
    % title('Voltage Profile','fontsize',12);
    set(gcf,'color','w');
    resetFigure(gcf);
    % ylim([0.8 1.2]);

    % data cursor mode
    dcm = datacursormode(fh);
    datacursormode on;
    set(dcm,'updatefcn',@dispBusName);

    o.Text.Command = 'Plot Profile Phases=ALL';
end
    function output_txt = dispBusName(obj,event_obj)
        % Display the position of the data cursor
        % obj          Currently not used (empty)
        % event_obj    Handle to event object
        % output_txt   Data cursor text string (string or cell array of strings).
        
        idx = find(dat.nodeDist == event_obj.Position(1));
        idy = find(dat.nodeV == event_obj.Position(2));
        id = intersect(idx,idy);
        if length(id) > 1
            dispBus = ['Buses: ' dat.nodeName{id(1)} ];
            for i = 2:length(id)
                dispBus = [dispBus '; ' dat.nodeName{id(i)}];
            end 
        else
            dispBus = ['Bus: ' dat.nodeName{id}];
        end
        output_txt = {dispBus, ['Dist: ' num2str(event_obj.Position(1))],...
            ['Volt: ' num2str(event_obj.Position(2))]};
    end
end