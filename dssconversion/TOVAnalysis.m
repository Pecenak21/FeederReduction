function TOVAnalysis(cc,p)

% analyses the SC behaviour of the provided feeders with and without
% presenece of PV
% INPUTS:
% cc - OpenDSS feeder structure
% p1 - feeder path. 

%% define the directory
[subdir subdir] = fileparts(fileparts(p));

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

%% Add monitors
n = dssmonitor('Name','Sub');% high voltage side
n.Element = ['Transformer.' cc.transformer(end).Name];
n.Mode = 0;
% Element = 1;
cc.monitor = n;

for yy=1:length(cc.line)
    n = dssmonitor('Name',[cc.line(yy).Name '.' cc.line(yy).bus1]);
    n.Element = ['Line.' cc.line(yy).Name];
    n.Mode = 0;
    n.Terminal = 1;
    cc.monitor(end+1) = n;
end

% SAVE MONITORS - OPTIONAL if you want to save all monitors as csv files
% after performing daily simulaitons. The DailySimulaiton.m runs the
% SCE_Centaur.dss file to export all monitors.
monitor_dir = [p(1:strfind(p,['/' subdir '.dss'])) 'Monitors'];
if ~exist([p(1:strfind(p,['/' subdir '.dss'])) 'Monitors'], 'dir');
    monitor_dir = [p(1:strfind(p,['/' subdir '.dss'])) 'Monitors'];
    mkdir(monitor_dir);
end

dfn = 'mon.dss';
for BB=2:length(cc.monitor)
    fid = fopen([monitor_dir '/' dfn], 'at');
    if(fid==-1), error('dsswrite:openfailed','Failed to open output file %s for writing!\nRemember to close open files before overwriting.',[savepath '/' dfn]); end
    str = ['export monitor ' cc.monitor(BB).Name];
    try
        fprintf(fid,'%s\n', str);
        fclose(fid);
    catch err
        warning('dsswrite:openfiles','Remember to close files before overwriting them!');
        rethrow(err);
    end
end

%% Save circuit
p = dsswrite(cc,subdir,1,subdir);

%% Run fault with PV (Centaur)
% [d_PV, vA_PV, vB_PV, vC_PV, vN_PV, iA_PV, iB_PV, iC_PV, iN_PV] = RunTOV(p,'PS0081_03265.1');%fault at the end of the feeder
% [d_PV, vA_PV, vB_PV, vC_PV, vN_PV, iA_PV, iB_PV, iC_PV, iN_PV] = RunTOV(p,'03265.1');%fault at the beginning of the feeder
% [d_PV, vA_PV, vB_PV, vC_PV, vN_PV, iA_PV, iB_PV, iC_PV, iN_PV] = RunTOV(p,'4120743E_03265_1.1');%fault at the mid-section of the feeder

%% Run fault with PV (Durox)
% [d_PV, vA_PV, vB_PV, vC_PV, vN_PV, iA_PV, iB_PV, iC_PV, iN_PV] = RunTOV(p,'ND63814968_05465.1');%fault at the end of the feeder
% [d_PV, vA_PV, vB_PV, vC_PV, vN_PV, iA_PV, iB_PV, iC_PV, iN_PV] = RunTOV(p,'05465.1');%fault at the beginning of the feeder
[d_PV, vA_PV, vB_PV, vC_PV, vN_PV, iA_PV, iB_PV, iC_PV, iN_PV] = RunTOV(p,'4258E_05465.1');%fault at the mid-section of the feeder

%% disable all pv units - by disabling the pv transformers, the pvs lose
%% their connection to the feeder. The way OpenDSS simulates this case is
%% by assigning the PV buses a NaN-Voltage for the fault condition. With the
%% transformers gone, it looks like the system has never had pv at all.
for i=1:length(cc.generator)% remove the transformers connecting the pvs to the feeder
    cc.transformer(i).enabled = 'no';
end
p = dsswrite(cc,subdir,1,subdir);

%% Run fault without PV (Centaur)
% [d_NoPV, vA_NoPV, vB_NoPV, vC_NoPV, vN_NoPV, iA_NoPV, iB_NoPV, iC_NoPV, iN_NoPV] = RunTOV(p,'PS0081_03265.1');%fault at the end of the feeder
% [d_NoPV, vA_NoPV, vB_NoPV, vC_NoPV, vN_NoPV, iA_NoPV, iB_NoPV, iC_NoPV, iN_NoPV] = RunTOV(p,'03265.1');%fault at the beginning of the feeder
% [d_NoPV, vA_NoPV, vB_NoPV, vC_NoPV, vN_NoPV, iA_NoPV, iB_NoPV, iC_NoPV, iN_NoPV] = RunTOV(p,'4120743E_03265_1.1');%fault at the mid-section of the feeder

%% Run fault without PV (Durox)
% [d_NoPV, vA_NoPV, vB_NoPV, vC_NoPV, vN_NoPV, iA_NoPV, iB_NoPV, iC_NoPV, iN_NoPV] = RunTOV(p,'ND63814968_05465.1');%fault at the end of the feeder
% [d_NoPV, vA_NoPV, vB_NoPV, vC_NoPV, vN_NoPV, iA_NoPV, iB_NoPV, iC_NoPV, iN_NoPV] = RunTOV(p,'05465.1');%fault at the beginning of the feeder
[d_NoPV, vA_NoPV, vB_NoPV, vC_NoPV, vN_NoPV, iA_NoPV, iB_NoPV, iC_NoPV, iN_NoPV] = RunTOV(p,'4258E_05465.1');%fault at the mid-section of the feeder

% % %% remove the results at the 120V level - they make the figure confusing
% % PvBusNames = cc.generator(:).bus1;
% % Names = regexprep(Names,'(\.\d+)+$','');
% % [aa bb] = ismember(Names, PvBusNames);
% % dA_PV=dA_PV(~aa); dB_PV=dB_PV(~aa); dC_PV=dC_PV(~aa); 
% % vA_PV=vA_PV(~aa); vB_PV=vB_PV(~aa); vC_PV=vC_PV(~aa);
% % % dN_PV=dN_PV(~aa); vN_PV=vN_PV(~aa);
% % vA_NoPV=vA_NoPV(~aa); vB_NoPV=vB_NoPV(~aa); vC_NoPV=vC_NoPV(~aa); 
% % % vN_NoPV=vN_NoPV(~aa);

%% Plot results for fault on phase A
plotPhases(vA_PV, vA_NoPV, d_PV, 'Voltage, pu', 'Phase A Voltage (Fault on Phase A)');%don't plot substation primary
% plotPhases(iA_PV, iA_NoPV, d_PV, 'Current, A', 'Phase A Current (Fault on Phase A)');
plotPhases(vB_PV, vB_NoPV, d_PV, 'Voltage, pu',  'Phase B Voltage (Fault on Phase A)');
% plotPhases(iB_PV, iB_NoPV, d_PV, 'Current, A', 'Phase B Current (Fault on Phase A)');
plotPhases(vC_PV, vC_NoPV, d_PV, 'Voltage, pu',  'Phase C Voltage (Fault on Phase A)');
% plotPhases(iC_PV, iC_NoPV, d_PV, 'Current, A', 'Phase C Current (Fault on Phase A)');
% plotPhases(vN_PV, vN_NoPV, d_PV, 'Voltage, pu',  'Neutral Voltage (Fault on Phase A)');
% plotPhases(iN_PV, iN_NoPV, d_PV, 'Current, A', 'Neutral Current (Fault on Phase A)');

%% RunTOV function
    function [dist, VA, VB, VC, VN, CurA, CurB, CurC, CurN] = RunTOV(p,FaultBus)
        o = actxserver('OpendssEngine.dss');
        c = o.ActiveCircuit;
        t = o.Text;
        m = c.Monitors;
        t.Command = 'Clear'; %Clear all circuits currently in memory.
        t.Command = ['Compile "' p '"']; %Reads the designated file name containing DSS commands and processes them as if they were entered directly into the command line.
        t.Command = 'solve mode=snap';
% % %         if exist('FaultBus', 'var')%if a fault bus is provided - execute the fault
% % % %             t.Command = 'solve mode=dynamics number=1 stepsize=0.00002';
% % %             t.Command = ['New Fault.f bus1=' FaultBus ' phases=1 r=0'];
% % % %             t.Command = ['New Fault.f bus1=PS0081_03265.1 phases=1 r=0'];
% % % %             t.Command = 'solve mode=snap';
% % %             t.Command = 'solve mode=direct';
% % %         end
            %% extract results
        t.Command = 'Sample';
%         DistA = o.ActiveCircuit.AllNodeDistancesByPhase(1)*(3280.8399)/1000;DistB = o.ActiveCircuit.AllNodeDistancesByPhase(2)*(3280.8399)/1000;DistC = o.ActiveCircuit.AllNodeDistancesByPhase(3)*(3280.8399)/1000;DistN = o.ActiveCircuit.AllNodeDistancesByPhase(4)*(3280.8399)/1000;
%         VoltA = o.ActiveCircuit.AllNodeVmagPUByPhase(1);VoltB = o.ActiveCircuit.AllNodeVmagPUByPhase(2);VoltC = o.ActiveCircuit.AllNodeVmagPUByPhase(3);VoltN = o.ActiveCircuit.AllNodeVmagPUByPhase(4);
%         Names = o.ActiveCircuit.AllNodeNamesByPhase(1);
        Dist = c.AllNodeDistances*(3280.8399)/1000;
        NodeNames = c.AllNodeNames';
        
        % OPTIONAL: exports all monitors to hard drive. See save monitors
        MonitorDirect = [p(1:strfind(p,['/' subdir '.dss'])) 'Monitors'];
        t.Command = ['cd ' MonitorDirect];
        t.Command = 'redirect mon.dss';
        
        CurA=zeros(1,461);  CurB=zeros(1,461); CurC=zeros(1,461); CurN=zeros(1,461);
        dist=zeros(1,461);
        VA=zeros(1,461); VB=zeros(1,461); VC=zeros(1,461); VN=zeros(1,461);
        i_A=1; i_B=1; i_C=1; i_N=1; i_T=1;
        Mon = m.First;
        while Mon>0,
            a = readMonitor(m.ByteStream);
            data = a.data;
            MonName = m.Name;
            k = regexp(MonName, '\.', 'split');
            if length(k)==1
                Mon = m.next;
                continue
            end
            bus = k{2};
            n = cell2mat(k([3:end])); 
            switch n
                case '123'
                    va=data(:,3);vb=data(:,5);vc=data(:,7);
                    ia=data(:,9);ia_ang=data(:,10);ib=data(:,11);ib_ang=data(:,12);ic=data(:,13);ic_ang=data(:,14);
                    CurA(:,i_A)=ia.*cosd(ia_ang); CurB(:,i_B) = ib.*cosd(ib_ang+120);CurC(:,i_C) = ic.*cosd(ic_ang-120);
                    bus = [bus '.1'];
                    [TF,LOC] = ismember(bus, NodeNames);
                    dist(:,i_T) = Dist(LOC);
                    VA(:,i_A) = va./(12000/sqrt(3)); VB(:,i_B) = vb./(12000/sqrt(3)); VC(:,i_C) = vc./(12000/sqrt(3));
                    i_A=i_A+1; i_B=i_B+1; i_C=i_C+1; i_T=i_T+1;
                 case '1234'%this case is only useful if neutrals are included
                    va=data(:,3);vb=data(:,5);vc=data(:,7);vn=data(:,9);
                    ia=data(:,11);ia_ang=data(:,12);ib=data(:,13);ib_ang=data(:,14);ic=data(:,15);ic_ang=data(:,16);in=data(:,17);in_ang=data(:,18);
                    CurA(:,i_A)=ia.*cosd(ia_ang); CurB(:,i_B) = ib.*cosd(ib_ang+120);CurC(:,i_C) = ic.*cosd(ic_ang-120);CurN(:,i_N) = in.*cosd(in_ang-120);
                    bus = [bus '.1'];
                    [TF,LOC] = ismember(bus, NodeNames);
                    dist(:,i_T) = Dist(LOC);
                    VA(:,i_A) = va./(12000/sqrt(3)); VB(:,i_B) = vb./(12000/sqrt(3)); VC(:,i_C) = vc./(12000/sqrt(3)); VN(:,i_N) = vn./(12000/sqrt(3));
                    i_A=i_A+1; i_B=i_B+1; i_C=i_C+1; i_N=i_N+1; i_T=i_T+1;
            end
            Mon = m.next;
        end        
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
