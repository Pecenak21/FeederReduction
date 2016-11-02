function FaultAnalysis(cc,FaultBus,FaultResistance)

% performs a single line to ground fault
% INPUTS:
% cc - OpenDSS feeder structure
% FaultBus - bus to host the fault (e.g. 03625.1)
% FaultResistance - fault resistance in ohms 

%% define the directory
savepath = pwd;
filename = cc.circuit.Name;

%% replace pv transformers with generators
len = length(cc.pvsystem);
for i=1:len
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
end
cc = rmfield(cc, 'pvsystem');

%% Save circuit
p = dsswrite(cc,filename,1,savepath);

%% Run fault with PV
[dA_PV, dB_PV, dC_PV, dN_PV, vA_PV, vB_PV, vC_PV, vN_PV, NamesA, NamesB, NamesC] = RunTOV(p,FaultBus,FaultResistance);

%% disable all pv units - by disabling the pv transformers, the pvs lose
%% their connection to the feeder. The way OpenDSS simulates this case is
%% by assigning the PV buses a NaN-Voltage for the fault condition. With the
%% transformers gone, it looks like the system has never had pv at all.
for i=1:length(cc.generator)% remove the transformers connecting the pvs to the feeder
    cc.transformer(i).enabled = 'no';
end
p = dsswrite(cc,filename,1,savepath);

%% Run fault without PV
[dA_NoPV, dB_NoPV, dC_NoPV, dN_NoPV, vA_NoPV, vB_NoPV, vC_NoPV, vN_NoPV] = RunTOV(p,FaultBus,FaultResistance);

%% remove the results at the 120V level - they make the figure confusing
PvBusNames = cc.generator(:).bus1;
Names = regexprep(NamesA,'(\.\d+)+$','');
[aa bb] = ismember(Names, PvBusNames);
dA_PV=dA_PV(~aa); vA_PV=vA_PV(~aa); vA_NoPV=vA_NoPV(~aa);
Names = regexprep(NamesB,'(\.\d+)+$','');
[aa bb] = ismember(Names, PvBusNames);
dB_PV=dB_PV(~aa); vB_PV=vB_PV(~aa); vB_NoPV=vB_NoPV(~aa);
Names = regexprep(NamesC,'(\.\d+)+$','');
[aa bb] = ismember(Names, PvBusNames);
dC_PV=dC_PV(~aa); vC_PV=vC_PV(~aa); vC_NoPV=vC_NoPV(~aa); 

%% Plot results
plotPhases(vA_PV(2:end), vA_NoPV(2:end), dA_PV(2:end), 'Voltage, pu', 'Phase A Voltage');%don't plot substation primary
plotPhases(vB_PV(2:end), vB_NoPV(2:end), dB_PV(2:end), 'Voltage, pu',  'Phase B Voltage');
plotPhases(vC_PV(2:end), vC_NoPV(2:end), dC_PV(2:end), 'Voltage, pu',  'Phase C Voltage');
% plotPhases(vN_PV(2:end), vN_NoPV(2:end), dN_PV(2:end), 'Voltage, pu',  'Neutral Voltage');

%% RunTOV function
    function [DistA, DistB, DistC, DistN, VoltA, VoltB, VoltC, VoltN, NamesA, NamesB, NamesC] = RunTOV(p,FaultBus, FaultResistance)
        o = actxserver('OpendssEngine.dss');
        t = o.Text;
        t.Command = 'Clear'; %Clear all circuits currently in memory.
        t.Command = ['Compile "' p '"']; %Reads the designated file name containing DSS commands and processes them as if they were entered directly into the command line.
        t.Command = 'solve mode=snap';
        if exist('FaultBus', 'var')%if a fault bus is provided - execute the fault
            t.Command = ['New Fault.f bus1=' FaultBus ' phases=1 r=' FaultResistance];
%             t.Command = 'solve mode=snap';
            t.Command = 'solve mode=direct';
        end
        %% extract results
        t.Command = 'Sample';
        DistA = o.ActiveCircuit.AllNodeDistancesByPhase(1)*(3280.8399)/1000;DistB = o.ActiveCircuit.AllNodeDistancesByPhase(2)*(3280.8399)/1000;DistC = o.ActiveCircuit.AllNodeDistancesByPhase(3)*(3280.8399)/1000;DistN = o.ActiveCircuit.AllNodeDistancesByPhase(4)*(3280.8399)/1000;
        VoltA = o.ActiveCircuit.AllNodeVmagPUByPhase(1);VoltB = o.ActiveCircuit.AllNodeVmagPUByPhase(2);VoltC = o.ActiveCircuit.AllNodeVmagPUByPhase(3);VoltN = o.ActiveCircuit.AllNodeVmagPUByPhase(4);
        NamesA = o.ActiveCircuit.AllNodeNamesByPhase(1);NamesB = o.ActiveCircuit.AllNodeNamesByPhase(2);NamesC = o.ActiveCircuit.AllNodeNamesByPhase(3);      
        o.delete;%close actserver
    end
%% plot function
    function [] = plotPhases(Fault, NoFault, Dist, ylabl, Title)
        dat_ = [Dist' Fault' Dist' NoFault' nan(length(Dist),2)];
        dat_ = reshape(dat_',2,numel(dat_)/2)';
        figure;
        h = plot(dat_(:,1), dat_(:,2),'r',Dist,Fault,'x',Dist,NoFault,'.');
        legend(h(2:3),{'PV','No PV'});
        xlabel('Distance, kft');
        ylabel(ylabl);
        title(Title);
        set(gcf,'nextplot','new');
        grid on;        
    end

end
