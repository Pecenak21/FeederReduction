function SCAnalysis(cc,p)

% analyses the SC behaviour of the provided feeders with and without
% presenece of PV
% INPUTS:
% c1 - OpenDSS feeder structure
% p1 - feeder path. 

%% define the directory
[subdir subdir] = fileparts(fileparts(p));

%% replace pv transformers with generators
len = length(cc.pvsystem);
for i=1:len
% %     bus1 = cc.transformer(i).Buses(1);
% %     bus2 = cc.transformer(i).Buses(2);
% %     name = cc.transformer(i).Name;
% %     cc.transformer(i).enabled = 'no';%remove transformer replaced by line
% %     n = dssline('Name',name);
% %     n.bus1 = bus1;
% %     n.bus2 = bus2;
% %     n.LineCode = 'BUS';
% %     n.Units = 'ft';
% %     n.Phases = 3;
% %     n.Length = 0.1;
% %     cc.line(end+1) = n;
    %%replace pv with generators
    name = cc.pvsystem(i).Name;
    m = dssgenerator('Name', name);
    m.Phases = 3;
    m.bus1 = cc.pvsystem(i).bus1;
    m.kv = 0.208;
    m.kw = 500;
    m.kvar = 0;
    m.Model = 7;
    m.Conn = 'delta';
    m.Vminpu = 0.01;
    m.Xdp = 0.069;
    m.balanced ='yes';
    if ~isfield(cc, 'generator')
        cc.generator = deal(m);
    else
        cc.generator(end+1) = m;
    end
    cc.pvsystem(i).Pmpp = 500;
    cc.pvsystem(i).kVA = 500;
    cc.pvsystem(i).Vminpu = 0.01;
    cc.pvsystem(i).LimitCurrent = 'yes';
end
cc = rmfield(cc, 'pvsystem');%remove pv object from the structure
p = dsswrite(cc,subdir,1,subdir);

%% load the fault study with pv
[y_PV, DistPV] = RunFault(p, '_fault_PV.csv');

%% extract results
ThreePhase_PV = cell2mat({y_PV.I3Phase}');%three phase fault currents
SinglePhase_PV = cell2mat({y_PV.I1Phase}');%three phase fault currents
LineLine_PV = cell2mat({y_PV.LL}');%three phase fault currents
BusNames = {y_PV.bus}';
PvBusNames = cc.generator(:).bus1;

%% remove the results at the 120V level - they make the figure confusing
[aa bb] = ismember(BusNames, PvBusNames);
ThreePhase_PV = ThreePhase_PV(~aa);
SinglePhase_PV = SinglePhase_PV(~aa);
LineLine_PV = LineLine_PV(~aa);
DistPV = DistPV(~aa);

%% disable all pv units - by disabling the pv transformers, the pvs lose
%% their connection to the feeder. The way OpenDSS simulates this case is
%% by assigning the PV buses a NaN-Voltage for the fault condition. With the
%% transformers gone, it looks like the system has never had pv at all.
for i=1:length(cc.generator)% remove the transformers connecting the pvs to the feeder
    cc.transformer(i).enabled = 'no';
%     cc.generator(i).enabled = 'no';
end
p = dsswrite(cc,subdir,1,subdir);

%% load the fault study without pv
[y_NoPV, Dist] = RunFault(p, '_fault_NoPV.csv');

%% extract results
ThreePhase_NoPV = cell2mat({y_NoPV.I3Phase}');%three phase fault currents
SinglePhase_NoPV = cell2mat({y_NoPV.I1Phase}');%three phase fault currents
LineLine_NoPV = cell2mat({y_NoPV.LL}');%three phase fault currents

%% remove the results at the 120V level
ThreePhase_NoPV = ThreePhase_NoPV(~aa);
SinglePhase_NoPV = SinglePhase_NoPV(~aa);
LineLine_NoPV = LineLine_NoPV(~aa);

%% plot results
plotFault(ThreePhase_PV([2:end]), ThreePhase_NoPV([2:end]), DistPV([2:end]), 'Three Phase Fault');
plotFault(SinglePhase_PV([2:end]), SinglePhase_NoPV([2:end]), DistPV([2:end]), 'Single Phase Fault');
plotFault(LineLine_PV([2:end]), LineLine_NoPV([2:end]), DistPV([2:end]), 'Line to Line Fault')

%% RunFault Function
    function [y, Dist] = RunFault(p, Result)
        o = actxserver('OpendssEngine.dss');
        t = o.Text;
        t.Command = 'Clear'; %Clear all circuits currently in memory.
        t.Command = ['Compile "' p '"']; %Reads the designated file name containing DSS commands and processes them as if they were entered directly into the command line.
        Dist = o.ActiveCircuit.AllBusDistances*(3280.8399)/1000; % conversion from km to kft
        t.Command = 'get editor'; olded = t.Result; if(isempty(olded)), olded = 'notepad.exe'; end
        t.Command = 'set editor="where.exe"'; % silence output
        t.Command = ['cd "' fileparts(p) '"/'];
        t.Command = 'Set mode=faultstudy';
        t.Command = 'Solve';
        t.Command = ['Export faultstudy ' lower(cc.circuit.Name) Result];
        y = faultread(t.Result);
        t.Command = 'show faults';
        t.Command = 'show currents';
        t.Command = ['cd "' pwd '"'];
        t.Command = ['set editor="' olded '"']; % unsilence output
        z = faultread([subdir '/' cc.circuit.Name '_FaultStudy.Txt']);
        o.delete;%close actserver
    end

%% plotFault Function
    function [] = plotFault(PV, NoPV, Dist, Title)
        dat_ = [Dist' PV Dist' NoPV nan(length(Dist),2)];
        dat_ = reshape(dat_',2,numel(dat_)/2)';
        figure;
        h = plot(dat_(:,1), dat_(:,2),'r',Dist,PV,'x',Dist,NoPV,'.');
        legend(h(2:3),{'PV','No PV'});
        xlabel('Distance, kft');
        ylabel('Current, A');
        title(Title);
        set(gcf,'nextplot','new');
        grid on;
    end
end
