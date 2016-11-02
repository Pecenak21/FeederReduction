function FaultAnalysis(cc,FaultBus,FaultResistance)

% performs a single line to ground fault
% INPUTS:
% cc - OpenDSS feeder structure
% FaultBus - bus to host the fault (e.g. 03625.1)
% FaultResistance - fault resistance in ohms 
cc=load(cc);
cc=cc.c1;

% [DSSStartOK, DSSObj, DSSText] = DSSStartup;
% 
% if DSSStartOK
%     DSSText.command='Compile (SCE_Centaur.dss)';
%     % Set up the interface variables
%     DSSCircuit=DSSObj.ActiveCircuit;
%     DSSSolution=DSSCircuit.Solution;
% end

%% define the directory
savepath = pwd;
filename = cc.circuit.Name;

%% replace pv transformers with generators
len = length(cc.pvsystem);
for i=1:len
%     cc.pvsystem(i).enabled='no';
    %%replace pv with generators
    name = cc.pvsystem(i).Name;
    m = dssgenerator('Name', name);
    m.Phases = 3;
    m.bus1 = cc.pvsystem(i).bus1;
    m.kv = 0.208;
    m.kw = 500;
    m.kvar = 0;
    m.Model = 7;
    m.Conn = 'wye';
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
p = dsswrite(cc,filename,1,savepath);

% %% modifying the system by changing line '90832664_03265_GS1336_3_03265'
% % i_line=GetIndex_struct(cc.line,'Name','90832664_03265_GS1336_3_03265');
% i_line=GetIndex_struct(cc.line,'Name','ND90832562_03265_GS1479_2_03265');
% % cc.line(i_line).enabled='no';
% cc.line(i_line).length=9999;
% % cc.line(i_line).R1=1;
% p = dsswrite(cc,filename,1,savepath);

% %% removing caps
% cc = rmfield(cc, 'capacitor');
% cc = rmfield(cc, 'capcontrol');
% p = dsswrite(cc,filename,1,savepath);

% %% removing all single phase load
% i_load=GetIndex_struct(cc.load,'Phases',1);
% cc.load(i_load)=[];
% p = dsswrite(cc,filename,1,savepath);
% 
% %% removing all three phase load
% i_load=GetIndex_struct(cc.load,'Phases',3);
% cc.load(i_load)=[];
% p = dsswrite(cc,filename,1,savepath);

% %% removing all loads
% cc = rmfield(cc, 'load');
% p = dsswrite(cc,filename,1,savepath);

% %% connection vsource to substation transformer (probably not necessary as default name is 'SourceBus')
% cc.circuit.bus1='SourceBus';
% p = dsswrite(cc,filename,1,savepath);

%% Run fault with PV
[dA_PV, dB_PV, dC_PV, dN_PV, vA_PV, vB_PV, vC_PV, vN_PV, NamesA, NamesB, NamesC] = RunTOV(p,FaultBus,FaultResistance);

%% disable all pv units - by disabling the pv transformers, the pvs lose their connection to the feeder. The way OpenDSS simulates this case is
%% by assigning the PV buses a NaN-Voltage for the fault condition. With the transformers gone, it looks like the system has never had pv at all.
% for i=1:length(cc.generator)% remove the transformers connecting the pvs to the feeder
%     cc.transformer(i).enabled = 'no';
% end
for i=1:length(cc.generator)% remove the transformers connecting the pvs to the feeder
    cc.transformer(i).enabled = 'no';
end
p = dsswrite(cc,filename,1,savepath);

%% Run fault without PV
[dA_NoPV, dB_NoPV, dC_NoPV, dN_NoPV, vA_NoPV, vB_NoPV, vC_NoPV, vN_NoPV] = RunTOV(p,FaultBus,FaultResistance);

%% remove the results at the 120V level - they make the figure confusing
PvBusNames = cc.generator(:).bus1;
% PvBusNames = cc.pvsystem(:).bus1;
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
close all;
plotPhases(vA_PV(2:end), vA_NoPV(2:end), dA_PV(2:end), 'Voltage, pu', 'A');%don't plot substation primary
plotPhases(vB_PV(2:end), vB_NoPV(2:end), dB_PV(2:end), 'Voltage, pu',  'B');
plotPhases(vC_PV(2:end), vC_NoPV(2:end), dC_PV(2:end), 'Voltage, pu',  'C');
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
        h = plot(dat_(:,1), dat_(:,2),'r',Dist,Fault,'^',Dist,NoFault,'o');
        set(gca,'FontSize',20);
        legend(h(2:3),{['PV, ' '\Phi' Title],['No PV, ' '\Phi' Title]},'FontSize',24);
        xlabel('Distance, kft','FontSize',26);
        ylabel(ylabl,'FontSize',26);
%         title(Title);
        set(h(2),'MarkerSize',5)
        set(h(3),'MarkerSize',5)
        set(h(2),'LineWidth',1)
        set(h(3),'LineWidth',1)
        set(h(2),'MarkerEdgeColor','r');
        set(h(2),'MarkerFaceColor','r');
        set(h(3),'MarkerEdgeColor','b');
        set(h(3),'MarkerFaceColor','b');
        if strcmp(Title,'B') || strcmp(Title,'C')
            ylim([1,1.05]);
        end
        saveas(gcf,[pwd '\' Title  '.tif'])
        saveas(gcf,[pwd '\' Title '.fig'])
        set(gcf,'nextplot','new');
        grid on;        
    end

end
