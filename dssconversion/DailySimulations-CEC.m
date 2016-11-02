function [] = DailySimulations(Commands, SimName, stepsize, SimPer, p, PVSys, NumPts)

% performes multiple simulaitons based on chosen stepsize and number of
% points.
% INPUTS: 
% Commands: any additional changes to the circuit on the fly (e.g. settign values, etc.)
% SimName: simulation title, will appear on the figures and figure names when saved
% stepsize: time will advance based on stepsize (e.g. 1h or 15m) for each solve
% NumPts: number of 'solves' (e.g. 24 for a daily simulaiton with 1h stepsize)
% p: path to main.dss 
% PVSys: names of all PV loadshapes in the circuit
% 3D Plots parameters
% z_saturation: z_saturation influences color spectrum, the higher the range, the farther away from extreme color
% z_range: range of values displayed, also determines z-axis numbering

disp('-- Running Daily Simulations --');

result_dir = [p(1:strfind(p,'/SCE_Centaur.dss')) 'results ' SimName];
step = str2double(stepsize(1:end-1));

if ~exist(result_dir, 'dir');
    result_dir = [p(1:strfind(p,'/SCE_Centaur.dss')) 'results ' SimName];
    mkdir(result_dir);
end

dssObj = actxserver('OpendssEngine.dss');
dssStart = dssObj.Start(0);
dssText = dssObj.Text;

%% OpendDSS simulaiton with PV included

dssText.Command = 'Clear'; %Clear all circuits currently in memory.
dssText.Command = ['Compile "' p '"']; %Reads the designated file name containing DSS commands and processes them as if they were entered directly into the command line. 
disp(['Simulating ' SimName ' with PV']);

if ~isempty(Commands) %any additional changes to the circuit on the fly (e.g. settign values, etc.)
    for x=1:1:length(Commands)
        dssText.Command = Commands{x};
    end
end

dssCircuit = dssObj.ActiveCircuit;
dssSolution = dssCircuit.Solution;

% dssEnergyMeters = dssCircuit.Meters;
dssMonitor = dssCircuit.Monitors;
% Dist = dssCircuit.AllNodeDistances*(3280.8399); % conversion from km to kft (not sure why OpenDSS presents it in km, when the line lengths are clearly identified to be in kft)
Dist = dssCircuit.AllNodeDistances; % Distance in km
NodeNames = dssCircuit.AllNodeNames';
% DistPhaseA = dssCircuit.AllNodeDistancesByPhase(1); DistPhaseB = dssCircuit.AllNodeDistancesByPhase(2); DistPhaseC = dssCircuit.AllNodeDistancesByPhase(3);
% NodeNamesPhaseA = dssCircuit.AllNodeNamesByPhase(1)'; NodeNamesPhaseB = dssCircuit.AllNodeNamesByPhase(2)'; NodeNamesPhaseC = dssCircuit.AllNodeNamesByPhase(3)'; 

dssSolution.MaxControlIterations=1000; %Maximum allowable control iterations (at 1000 max iterations, the simulaiton might take a while but the chances to get a solution increase)
dssSolution.maxiterations=1000; %Maximum allowable number of iterations for the circuit solution (not counting control iterations)

%{OFF | STATIC |EVENT | TIME} - STATIC = Time does not advance. Control actions are executed in order of shortest time to act until all actions are cleared from the control queue. 
%Use this mode for power flow solutions which may require several regulator tap changes per solution. This is the default for the standard Snapshot mode as well as 
%Daily and Yearly simulations where the stepsize is typically greater than 15 min.
dssText.Command = 'Set controlmode = static';
dssText.Command = ['Set mode = daily stepsize = ' stepsize];% time will advance based on stepsize (e.g. 1h or 15m) for each solve
dssText.Command = 'Set number = 1';% will stop after each solve

% This function initializes iteration counters, etc that occur at the beginning of the SolveSnap function. Invoke this method if you are implementing your own solution
% process. Calls the SnapShotInit function in the DSS.
dssSolution.InitSnap; 

% Present time expressed as a floating point number in units of hours.
dssSolution.dblHour = 0.0;

Volt_MaxMin_PV=zeros(NumPts,2); LossTotal_PV=zeros(NumPts,2); LossLine_PV=zeros(NumPts,2);
Volt_PV=zeros(NumPts,length(dssCircuit.AllBusVmagPu)); 
% TotalPower_PV=zeros(NumPts,2);

i=1;
%tic
while (dssSolution.dblHour < SimPer)
    dssSolution.Solve;
    if dssSolution.Converged
        Volt_MaxMin_PV(i,:) = [max(dssCircuit.AllBusVmagPu) min(dssCircuit.AllBusVmagPu)];
        LossTotal_PV(i,:) = dssCircuit.Losses/(1e+06);
        LossLine_PV(i,:) = dssCircuit.LineLosses/(1e+06);
        Volt_PV(i,:) = dssCircuit.AllBusVmagPu;
            %% Returns the total power in kW and kvar supplied to the circuit 
            %% by all Vsource and Isource objects. Does not include Generator 
            %% objects. Returned as a two-element array of doubles........
            %% Have to be carefull though - this numbers include the
            %% substation transformer losses as well. If no details about
            %% the substation transformer are known it is better to use the
            %% kw and kvars from the line monitor coming out of the
            %% substation transformer secondary. Otherwise, if bad
            %% transformer parameters were chosen these numbers can be
            %% quite missleading.
%         TotalPower_PV(i,:) = dssCircuit.TotalPower/(1e+03)*(-1);
        i=i+1;
    else
        disp(['System did not converge for ''with PV'' case (NICHT GUT) at hour: ' num2str(dssSolution.dblHour)]);
    end
    
end
%toc
% keyboard
EventLog = dssSolution.EventLog;
[Event_PV Reg_PV Cap_PV] = EventLogEvaluation(EventLog, i-1, step);

% read all the monitors in the circuit
disp('- Evaluating EventLog');
CurA_PV=zeros(NumPts,464);  CurB_PV=zeros(NumPts,464); CurC_PV=zeros(NumPts,461);
kWA_PV=zeros(NumPts,464); kWB_PV=zeros(NumPts,464); kWC_PV=zeros(NumPts,461); kW_PV=zeros(NumPts,464);
kVarA_PV=zeros(NumPts,464); kVarB_PV=zeros(NumPts,464); kVarC_PV=zeros(NumPts,461); kVar_PV=zeros(NumPts,464);
distA=zeros(1,464); distB=zeros(1,464); distC=zeros(1,461); dist=zeros(1,464);
VA_PV=zeros(NumPts,464); VB_PV=zeros(NumPts,464); VC_PV=zeros(NumPts,461);
i_A=1; i_B=1; i_C=1; i_T=1;
Mon = dssMonitor.First;
%tic
while Mon>0,
    a = readMonitor(dssMonitor.ByteStream);
    data = a.data;
    MonName = dssMonitor.Name;
    m = regexp(MonName, '\.', 'split');
    if length(m)==1
        Mon = dssMonitor.next;
        continue
    end
    bus = m{2};
    n = cell2mat(m([3:end]));
    switch n
        case '12'
            va=data(:,3);va_ang=data(:,4);vb=data(:,5);vb_ang=data(:,6);
            ia=data(:,7);ia_ang=data(:,8);ib=data(:,9);ib_ang=data(:,10);
            CurA_PV(:,i_A) = ia.*cosd(ia_ang);CurB_PV(:,i_B) = ib.*cosd(ib_ang+120);
            kWA_PV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;kWB_PV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;
            kW_PV(:,i_T) = kWA_PV(:,i_A)+kWB_PV(:,i_B);
            kVarA_PV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarB_PV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;
            kVar_PV(:,i_T) = kVarA_PV(:,i_A)+kVarB_PV(:,i_B);
            bus = [bus '.1'];
            [TF,LOC] = ismember(bus, NodeNames);
            distA(:,i_A) = Dist(LOC); distB(:,i_B) = Dist(LOC); dist(:,i_T) = Dist(LOC);
            VA_PV(:,i_A) = va./(12000/sqrt(3)); VB_PV(:,i_B) = vb./(12000/sqrt(3));
            i_A=i_A+1; i_B=i_B+1; i_T=i_T+1;
        case '13'
            va=data(:,3);va_ang=data(:,4);vc=data(:,5);vc_ang=data(:,6);
            ia=data(:,7);ia_ang=data(:,8);ic=data(:,9);ic_ang=data(:,10);
            CurA_PV(:,i_A) = ia.*cosd(ia_ang);CurC_PV(:,i_C) = ic.*cosd(ic_ang+120);
            kWA_PV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;kWC_PV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
            kW_PV(:,i_T) = kWA_PV(:,i_A)+kWC_PV(:,i_C);
            kVarA_PV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarC_PV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
            kVar_PV(:,i_T) = kVarA_PV(:,i_A)+kVarC_PV(:,i_C);
            bus = [bus '.1'];
            [TF,LOC] = ismember(bus, NodeNames);
            distA(:,i_A) = Dist(LOC); distC(:,i_C) = Dist(LOC);dist(:,i_T) = Dist(LOC);
            VA_PV(:,i_A) = va./(12000/sqrt(3)); VC_PV(:,i_C) = vc./(12000/sqrt(3));
            i_A=i_A+1; i_C=i_C+1; i_T=i_T+1;
        case '23'
            vb=data(:,3);vb_ang=data(:,4);vc=data(:,5);vc_ang=data(:,6);
            ib=data(:,7);ib_ang=data(:,8);ic=data(:,9);ic_ang=data(:,10);
            CurB_PV(:,i_B) = ib.*cosd(ib_ang);CurC_PV(:,i_C) = ic.*cosd(ic_ang+120);
            kWB_PV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;kWC_PV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
            kW_PV(:,i_T) = kWB_PV(:,i_B)+kWC_PV(:,i_C);
            kVarB_PV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;kVarC_PV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
            kVar_PV(:,i_T) = kVarB_PV(:,i_B)+kVarC_PV(:,i_C);
            bus = [bus '.2'];
            [TF,LOC] = ismember(bus, NodeNames);
            distB(:,i_B) = Dist(LOC); distC(:,i_C) = Dist(LOC); dist(:,i_T) = Dist(LOC);
            VB_PV(:,i_B) = vb./(12000/sqrt(3)); VC_PV(:,i_C) = vc./(12000/sqrt(3));
            i_B=i_B+1; i_C=i_C+1; i_T=i_T+1;
        case '123'
            va=data(:,3);va_ang=data(:,4);vb=data(:,5);vb_ang=data(:,6);vc=data(:,7);vc_ang=data(:,8);
            ia=data(:,9);ia_ang=data(:,10);ib=data(:,11);ib_ang=data(:,12);ic=data(:,13);ic_ang=data(:,14);
            CurA_PV(:,i_A)=ia.*cosd(ia_ang); CurB_PV(:,i_B) = ib.*cosd(ib_ang+120);CurC_PV(:,i_C) = ic.*cosd(ic_ang-120);
            kWA_PV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)./1000;kWB_PV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;kWC_PV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
            kW_PV(:,i_T) = kWA_PV(:,i_A)+kWB_PV(:,i_B)+kWC_PV(:,i_C);
            kVarA_PV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarB_PV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;kVarC_PV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
            kVar_PV(:,i_T) = kVarA_PV(:,i_A)+kVarB_PV(:,i_B)+kVarC_PV(:,i_C);
            bus = [bus '.1'];
            [TF,LOC] = ismember(bus, NodeNames);
            distA(:,i_A) = Dist(LOC); distB(:,i_B) = Dist(LOC); distC(:,i_C) = Dist(LOC);dist(:,i_T) = Dist(LOC);
            VA_PV(:,i_A) = va./(12000/sqrt(3)); VB_PV(:,i_B) = vb./(12000/sqrt(3)); VC_PV(:,i_C) = vc./(12000/sqrt(3));
            i_A=i_A+1; i_B=i_B+1; i_C=i_C+1; i_T=i_T+1;
        case '1'
            va=data(:,3);va_ang=data(:,4);
            ia=data(:,5);ia_ang=data(:,6);
            CurA_PV(:,i_A) = ia.*cosd(ia_ang);
            kWA_PV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;
            kW_PV(:,i_T) = kWA_PV(:,i_A);
            kVarA_PV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;
            kVar_PV(:,i_T) = kVarA_PV(:,i_A);
            bus = [bus '.1'];
            [TF,LOC] = ismember(bus, NodeNames);
            distA(:,i_A) = Dist(LOC); dist(:,i_T) = Dist(LOC);
            VA_PV(:,i_A) = va./(12000/sqrt(3));
            i_A=i_A+1;  i_T=i_T+1;
        case '2'
            vb=data(:,3);vb_ang=data(:,4);
            ib=data(:,5);ib_ang=data(:,6);
            CurB_PV(:,i_B) = ib.*cosd(ib_ang+120);
            kWB_PV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;
            kW_PV(:,i_T) = kWB_PV(:,i_B);
            kVarB_PV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;
            kVar_PV(:,i_T) = kVarB_PV(:,i_B);
            bus = [bus '.2'];
            [TF,LOC] = ismember(bus, NodeNames);
            distB(:,i_B) = Dist(LOC); dist(:,i_T) = Dist(LOC);
            VB_PV(:,i_B) = vb./(12000/sqrt(3));
            i_B=i_B+1;  i_T=i_T+1;
        case '3'
            vc=data(:,3);vc_ang=data(:,4);
            ic=data(:,5);ic_ang=data(:,6);
            CurC_PV(:,i_C) = ic.*cosd(ic_ang-120);
            kWC_PV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
            kW_PV(:,i_T) = kWC_PV(:,i_C);
            kVarC_PV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
            kVar_PV(:,i_T) = kVarC_PV(:,i_C);
            bus = [bus '.3'];
            [TF,LOC] = ismember(bus, NodeNames);
            distC(:,i_C) = Dist(LOC);dist(:,i_T) = Dist(LOC);
            VC_PV(:,i_C) = vc./(12000/sqrt(3));
            i_C=i_C+1;  i_T=i_T+1;
    end
    
    Mon = dssMonitor.next;
end
%toc
% OPTIONAL: exports all monitors to hard drive. See save monitors section
% in main_f520.m
% MonitorDirect = [p(1:strfind(p,'/SCE_Centaur.dss')) 'Monitors'];
% dssText.Command = ['cd ' MonitorDirect];
% dssText.Command = 'redirect Centaur_mon.dss';

%% OpendDSS simulaiton with No PV included

dssText.Command = 'Clear';% clears out the mamory from the privious simulation
dssText.Command = ['Compile "' p '"'];% compiles the existing circuit, anew
disp(['Simulating ' SimName ' without PV']);

if ~isempty(Commands)% it is important to have the same set of desired changes for this system as well
    for x=1:1:length(Commands)
        dssText.Command = Commands{x};
    end
end

for a=1:length(PVSys)
    dssText.Command = ['Loadshape.' PVSys{a} '.Mult = (' num2str(zeros(1,NumPts)) ')'];
% %     checks if the pv loadshapes were set to 0 - just to make sure
%     dssText.Command = ['? Loadshape.' PVSys{a} '.Mult'];
%     dssText.Result
end

dssCircuit = dssObj.ActiveCircuit;
dssSolution = dssCircuit.Solution;

dssSolution.MaxControlIterations=1000;
dssSolution.maxiterations=1000;

dssText.Command = 'Set controlmode = static';
dssText.Command = ['Set mode = daily stepsize = ' stepsize];
dssText.Command = 'Set number = 1';% will stop after each solve

dssSolution.InitSnap;
dssSolution.dblHour = 0.0;
 
Volt_MaxMin_NoPV=zeros(NumPts,2); LossTotal_NoPV=zeros(NumPts,2); LossLine_NoPV=zeros(NumPts,2);
Volt_NoPV=zeros(NumPts,length(dssCircuit.AllBusVmagPu)); 
% TotalPower_NoPV=zeros(NumPts,2);

i=1;
%tic
while (dssSolution.dblHour < SimPer)
    dssSolution.Solve;
    if dssSolution.Converged
        Volt_MaxMin_NoPV(i,:) = [max(dssCircuit.AllBusVmagPu) min(dssCircuit.AllBusVmagPu)];
        LossTotal_NoPV(i,:) = dssCircuit.Losses/(1e+06);
        Volt_NoPV(i,:) = dssCircuit.AllBusVmagPu;
        LossLine_NoPV(i,:) = dssCircuit.LineLosses/(1e+06);
            %% Returns the total power in kW and kvar supplied to the circuit 
            %% by all Vsource and Isource objects. Does not include Generator 
            %% objects. Returned as a two-element array of doubles........
            %% Have to be carefull though - this numbers include the
            %% substation transformer losses as well. If no details about
            %% the substation transformer are known it is better to use the
            %% kw and kvars from the line monitor coming out of the
            %% substation transformer secondary. Otherwise, if bad
            %% transformer parameters were chosen these numbers can be
            %% quite missleading.
%         TotalPower_NoPV(i,:) = dssCircuit.TotalPower/(1e+03)*(-1);
        i=i+1;
    else
        disp(['System did not converge for ''no PV'' case (NICHT GUT) at hour: ' num2str(dssSolution.dblHour)]);
    end
    
end
%toc
EventLog = dssSolution.EventLog;
[Event_NoPV Reg_NoPV Cap_NoPV] = EventLogEvaluation(EventLog, i-1, step);

disp('- Evaluating EventLog');
CurA_NoPV=zeros(NumPts,464);  CurB_NoPV=zeros(NumPts,464); CurC_NoPV=zeros(NumPts,461);
kWA_NoPV=zeros(NumPts,464); kWB_NoPV=zeros(NumPts,464); kWC_NoPV=zeros(NumPts,461); kW_NoPV=zeros(NumPts,464);
kVarA_NoPV=zeros(NumPts,464); kVarB_NoPV=zeros(NumPts,464); kVarC_NoPV=zeros(NumPts,461); kVar_NoPV=zeros(NumPts,464);
% distA=zeros(1,2032); distB=zeros(1,2047); distC=zeros(1,2020); dist=zeros(1,2454);
VA_NoPV=zeros(NumPts,464); VB_NoPV=zeros(NumPts,464); VC_NoPV=zeros(NumPts,461);
i_A=1; i_B=1; i_C=1; i_T=1;
Mon = dssMonitor.First;
%tic
while Mon>0,
    a = readMonitor(dssMonitor.ByteStream);
    data = a.data;
    MonName = dssMonitor.Name;
    m = regexp(MonName, '\.', 'split');
    if length(m)==1
        Mon = dssMonitor.next;
        continue
    end
    bus = m{2};
    n = cell2mat(m([3:end]));
    switch n
        case '12'
            va=data(:,3);va_ang=data(:,4);vb=data(:,5);vb_ang=data(:,6);
            ia=data(:,7);ia_ang=data(:,8);ib=data(:,9);ib_ang=data(:,10);
            CurA_NoPV(:,i_A) = ia.*cosd(ia_ang);CurB_NoPV(:,i_B) = ib.*cosd(ib_ang+120);
            kWA_NoPV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;kWB_NoPV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;
            kW_NoPV(:,i_T) = kWA_NoPV(:,i_A)+kWB_NoPV(:,i_B);
            kVarA_NoPV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarB_NoPV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;
            kVar_NoPV(:,i_T) = kVarA_NoPV(:,i_A)+kVarB_NoPV(:,i_B);
            bus = [bus '.1'];
            VA_NoPV(:,i_A) = va./(12000/sqrt(3));VB_NoPV(:,i_B) = vb./(12000/sqrt(3));
            i_A=i_A+1; i_B=i_B+1; i_T=i_T+1;
        case '13'
            va=data(:,3);va_ang=data(:,4);vc=data(:,5);vc_ang=data(:,6);
            ia=data(:,7);ia_ang=data(:,8);ic=data(:,9);ic_ang=data(:,10);
            CurA_NoPV(:,i_A) = ia.*cosd(ia_ang);CurC_NoPV(:,i_C) = ic.*cosd(ic_ang+120);
            kWA_NoPV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;kWC_NoPV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
            kW_NoPV(:,i_T) = kWA_NoPV(:,i_A)+kWC_NoPV(:,i_C);
            kVarA_NoPV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarC_NoPV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
            kVar_NoPV(:,i_T) = kVarA_NoPV(:,i_A)+kVarC_NoPV(:,i_C);
            bus = [bus '.1'];
            VA_NoPV(:,i_A) = va./(12000/sqrt(3));VC_NoPV(:,i_C) = vc./(12000/sqrt(3));
            i_A=i_A+1; i_C=i_C+1; i_T=i_T+1;
        case '23'
            vb=data(:,3);vb_ang=data(:,4);vc=data(:,5);vc_ang=data(:,6);
            ib=data(:,7);ib_ang=data(:,8);ic=data(:,9);ic_ang=data(:,10);
            CurB_NoPV(:,i_B) = ib.*cosd(ib_ang);CurC_NoPV(:,i_C) = ic.*cosd(ic_ang+120);
            kWB_NoPV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;kWC_NoPV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
            kW_NoPV(:,i_T) = kWB_NoPV(:,i_B)+kWC_NoPV(:,i_C);
            kVarB_NoPV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;kVarC_NoPV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
            kVar_NoPV(:,i_T) = kVarB_NoPV(:,i_B)+kVarC_NoPV(:,i_C);
            bus = [bus '.2'];
            VB_NoPV(:,i_B) = vb./(12000/sqrt(3));VC_NoPV(:,i_C) = vc./(12000/sqrt(3));
            i_B=i_B+1; i_C=i_C+1; i_T=i_T+1;
        case '123'
            va=data(:,3);va_ang=data(:,4);vb=data(:,5);vb_ang=data(:,6);vc=data(:,7);vc_ang=data(:,8);
            ia=data(:,9);ia_ang=data(:,10);ib=data(:,11);ib_ang=data(:,12);ic=data(:,13);ic_ang=data(:,14);
            CurA_NoPV(:,i_A) = ia.*cosd(ia_ang);CurB_NoPV(:,i_B) = ib.*cosd(ib_ang+120);CurC_NoPV(:,i_C) = ic.*cosd(ic_ang-120);
            kWA_NoPV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;kWB_NoPV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;kWC_NoPV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
            kW_NoPV(:,i_T) = kWA_NoPV(:,i_A)+kWB_NoPV(:,i_B)+kWC_NoPV(:,i_C);
            kVarA_NoPV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarB_NoPV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;kVarC_NoPV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
            kVar_NoPV(:,i_T) = kVarA_NoPV(:,i_A)+kVarB_NoPV(:,i_B)+kVarC_NoPV(:,i_C);
            bus = [bus '.1'];
            VA_NoPV(:,i_A) = va./(12000/sqrt(3));VB_NoPV(:,i_B) = vb./(12000/sqrt(3));VC_NoPV(:,i_C) = vc./(12000/sqrt(3));
            i_A=i_A+1; i_B=i_B+1; i_C=i_C+1; i_T=i_T+1;
        case '1'
            va=data(:,3);va_ang=data(:,4);
            ia=data(:,5);ia_ang=data(:,6);
            CurA_NoPV(:,i_A) = ia.*cosd(ia_ang);
            kWA_NoPV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;
            kW_NoPV(:,i_T) = kWA_NoPV(:,i_A);
            kVarA_NoPV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;
            kVar_NoPV(:,i_T) = kVarA_NoPV(:,i_A);
            bus = [bus '.1'];
            VA_NoPV(:,i_A) = va./(12000/sqrt(3));
            i_A=i_A+1; i_T=i_T+1;
        case '2'
            vb=data(:,3);vb_ang=data(:,4);
            ib=data(:,5);ib_ang=data(:,6);
            CurB_NoPV(:,i_B) = ib.*cosd(ib_ang+120);
            kWB_NoPV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;
            kW_NoPV(:,i_T) = kWB_NoPV(:,i_B);
            kVarB_NoPV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;
            kVar_NoPV(:,i_T) = kVarB_NoPV(:,i_B);
            bus = [bus '.2'];
            VB_NoPV(:,i_B) = vb./(12000/sqrt(3));
            i_B=i_B+1; i_T=i_T+1;
        case '3'
            vc=data(:,3);vc_ang=data(:,4);
            ic=data(:,5);ic_ang=data(:,6);
            CurC_NoPV(:,i_C) = ic.*cosd(ic_ang-120);
            kWC_NoPV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
            kW_NoPV(:,i_T) = kWC_NoPV(:,i_C);
            kVarC_NoPV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
            kVar_NoPV(:,i_T) = kVarC_NoPV(:,i_C);
            bus = [bus '.3'];
            VC_NoPV(:,i_C) = vc./(12000/sqrt(3));
            i_C=i_C+1; i_T=i_T+1;
    end
    
    Mon = dssMonitor.next;
end
%toc
% save([result_dir '/' SimName '_results.mat'], 'LossTotal_PV', 'LossTotal_NoPV', 'TotalPower_PV', 'TotalPower_NoPV');
% OPTIONAL: exports all monitors to hard drive. See save monitors section
% in main_f520.m
% MonitorDirect = [p(1:strfind(p,'/f520.dss')) 'Monitors'];
% dssText.Command = ['cd ' MonitorDirect];
% dssText.Command = 'redirect f520_mon.dss';

%% 2D Plots
% time = i-1;
% xrange = 1:(i-1);
% if NumPts > 2000
%     time=time/120;
%     xrange = xrange./120;
% end

StartDate = datenum('00:00:00','HH:MM:SS');
EndDate = datenum('23:59:45','HH:MM:SS');
xrange = (StartDate:(EndDate-StartDate)/(NumPts-1):EndDate);

disp(['Plotting results for ' SimName]);
%tic
% if the visibility for the 2D plots is set to 'off', Matlab does a poor
% job saving .fig files - they don't open for editing, anymore. the
% simulaiton gets slower but at least more editing can be done later on.

%%plot substation voltages. assume the substation secondary monitor is the first one in the list
f = figure; hold on;
h = plot(xrange, VA_NoPV(:,1), 'r', xrange, VB_NoPV(:,1), 'g', xrange, VC_NoPV(:,1), 'b');
legend(h,{'Phase A', 'Phase B','Phase C'},'fontsize',20);
title('Substation Phase Voltages - No PV','fontsize',20)
xlabel('Time of Day','fontsize',20)
ylabel('Voltages (pu)','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylim([0.9 1.1])
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

f = figure; hold on;
h = plot(xrange, VA_PV(:,1), 'r', xrange, VB_PV(:,1), 'g', xrange, VC_PV(:,1), 'b');
legend(h,{'Phase A', 'Phase B','Phase C'},'fontsize',20);
title('Substation Phase Voltages - With PV','fontsize',20)
xlabel('Time of Day','fontsize',20)
ylabel('Voltages (pu)','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylim([0.9 1.1])
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

f = figure; hold on;
h = plot(xrange, Volt_MaxMin_PV(:,1), 'r', xrange, Volt_MaxMin_NoPV(:,1),...
    'g', xrange, Volt_MaxMin_PV(:,2), 'b', xrange, Volt_MaxMin_NoPV(:,2), 'k');
legend(h,{'Max (PV)', 'Max (No PV)','Min (PV)', 'Min (No PV)'},'fontsize',20);
title('Maximum and Minimum Bus Voltages','fontsize',20)
xlabel('Time of Day','fontsize',20)
ylabel('Voltages (pu)','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylim([0.9 1.1])
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

f = figure; hold on;
h = plot(xrange, LossTotal_PV(:,1), 'r', xrange, LossTotal_NoPV(:,1),...
    'g', xrange, LossTotal_PV(:,2), 'b', xrange, LossTotal_NoPV(:,2), 'k');
legend(h,{'MW (PV)', 'MW (No PV)', 'MVar (PV)', 'MVar (No PV)'},'fontsize',20)
title('Total Losses','fontsize',20)
xlabel('Time of Day','fontsize',20)
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylim([0 0.5])
ylabel('Losses','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

f = figure; hold on;
h = plot(xrange, LossLine_PV(:,1), 'r', xrange, LossLine_NoPV(:,1),...
    'g', xrange, LossLine_PV(:,2), 'b', xrange, LossLine_NoPV(:,2), 'k');
legend(h,{'MW (PV)', 'MW (No PV)', 'MVar (PV)', 'MVar (No PV)'},'fontsize',20)
title('Total Line Losses','fontsize',20)
xlabel('Time of Day','fontsize',20)
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylabel('Losses','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

f = figure; hold on;
    %% Instead of using the total power data provided by OpenDSS use the kw
    %% and kvar data obtained from the monitor connected to the line coming
    %% out of the substation. See explanation above (OpenDSS similation loop)
% % % h = plot(xrange, TotalPower_PV(:,1), 'r', xrange, TotalPower_NoPV(:,1),...
% % %     'g', xrange, TotalPower_PV(:,2).*(-1), 'b', xrange, TotalPower_NoPV(:,2).*(-1), 'k');
h = plot(xrange, kW_PV(:,1)./1000, 'r', xrange, kW_NoPV(:,1)./1000,...
    'g', xrange, kVar_PV(:,1)./1000, 'b', xrange, kVar_NoPV(:,1)./1000, 'k');
legend(h,{'MW (PV)', 'MW (No PV)', 'MVar (PV)', 'MVar (No PV)'},'fontsize',20)
title('Substation Power','fontsize',20)
xlabel('Time of Day','fontsize',20)
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylim([-5 10])
ylabel('Power','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

f = figure; hold on;
h = plot(xrange, Reg_PV((2:end),1), 'r*-', xrange, Reg_PV((2:end),2), 'g*-');
legend(h,{['I (' num2str(Reg_PV(1,1)) ')'], ['Sub (' num2str(Reg_PV(1,2)) ')']})
title('Voltage Regulator Events (With PV)','fontsize',20)
xlabel('Time of Day','fontsize',20)
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylabel('# Events','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

f = figure; hold on;
h = plot(xrange, Reg_NoPV((2:end),1), 'r*-', xrange, Reg_NoPV((2:end),2), 'g');
legend(h,{['I (' num2str(Reg_NoPV(1,1)) ')'],['Sub (' num2str(Reg_NoPV(1,2)) ')']})
title('Voltage Regulator Events (No PV)','fontsize',20)
xlabel('Time of Day','fontsize',20)
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylabel('# Events','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

f = figure; hold on;
h = plot(xrange, Cap_PV((2:end),1), 'r');
legend(h,{['I (' num2str(Cap_PV(1,1)) ')']})
title('Capacitor Bank Events (With PV)','fontsize',20)
xlabel('Time of Day','fontsize',20)
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylabel('# Events','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

f = figure; hold on;
h = plot(xrange, Cap_NoPV((2:end),1), 'r');
legend(h,{['I (' num2str(Cap_NoPV(1,1)) ')']})
title('Capacitor Bank Events (No PV)','fontsize',20)
xlabel('Time of Day','fontsize',20)
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylabel('# Events','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

RegTotal_PV = sum(Reg_PV(2:end,end));
RegTotal_NoPV = sum(Reg_NoPV(2:end,end));
f = figure; hold on;
h = plot(xrange, Reg_PV(2:end,end), 'r', xrange, Reg_NoPV(2:end,end), 'g');
legend(h,{['Total_P_V (' num2str(RegTotal_PV) ')'], ['Total_N_o_P_V (' num2str(RegTotal_NoPV) ')']})
title('Voltage Regulator Events - Total','fontsize',20)
xlabel('Time of Day','fontsize',20)
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylabel('# Events','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

CapTotal_PV = sum(Cap_PV(2:end,end));
CapTotal_NoPV = sum(Cap_NoPV(2:end,end));
f = figure; hold on;
h = plot(xrange, Cap_PV(2:end,end), 'r', xrange, Cap_NoPV(2:end,end), 'g');
legend(h,{['Total_P_V (' num2str(CapTotal_PV) ')'], ['Total_N_o_P_V (' num2str(CapTotal_NoPV) ')']})
title('Capacitor Bank Events - Total','fontsize',20)
xlabel('Time of Day','fontsize',20)
datetick('x','HH:MM','keeplimits', 'keepticks')
xlim([min(xrange) max(xrange)])
ylabel('# Events','fontsize',20)
set(0,'DefaultAxesFontSize', 20)
grid on;
    %% set figure to fit the screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
    %% save figure
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

close all

%toc
%% 3D Plots
% 
Directory_Base_Output = result_dir;
TimeMin=1;
TimeMax=i-1;
Edge=0; % 0 -> no contour line (preferred for a large number of cases

X_Label_3D(1)={'Distance from'};
X_Label_3D(2)={'Substation (km)'};
Y_Label_3D(1)={'Time of Day'};
Z_Label_3D(1)={'Voltage, pu'};
X_Label_2D(1)={'Distance from Substation (km)'};
Y_Label_2D(1)={'Time of Day'};
Z_Label_2D(1)={'Voltage, pu'};

Label_FontSize_2D=9;
Label_FontSize_3D=9;

DTick = 'HH:MM';

z_saturation=[0.90 1.15];% z_saturation influences color spectrum, the higher the range, the farther away from extreme color
z_range=(0.9:0.05:1.10); % range of values displayed, also determines z-axis numbering
% z_saturation=[0.90 1.15];
% z_range=(0.5:0.05:1.10);

if NumPts > 2000 % downsize the number of points to be plotted
    nth = 10;% plot every nth point
    Volt_NoPV=Volt_NoPV(1:nth:end,:);
    Volt_PV=Volt_PV(1:nth:end,:);
    VA_NoPV=VA_NoPV(1:nth:end,:);
    VA_PV=VA_PV(1:nth:end,:);
    VB_NoPV=VB_NoPV(1:nth:end,:);
    VB_PV=VB_PV(1:nth:end,:);
    VC_NoPV=VC_NoPV(1:nth:end,:);
    VC_PV=VC_PV(1:nth:end,:);
    kWA_NoPV=kWA_NoPV(1:nth:end,:);
    kWA_PV=kWA_PV(1:nth:end,:);
    kWB_NoPV=kWB_NoPV(1:nth:end,:);
    kWB_PV=kWB_PV(1:nth:end,:);
    kWC_NoPV=kWC_NoPV(1:nth:end,:);
    kWC_PV=kWC_PV(1:nth:end,:);
    kW_NoPV=kW_NoPV(1:nth:end,:);
    kW_PV=kW_PV(1:nth:end,:);
    kVarA_NoPV=kVarA_NoPV(1:nth:end,:);
    kVarA_PV=kVarA_PV(1:nth:end,:);
    kVarB_NoPV=kVarB_NoPV(1:nth:end,:);
    kVarB_PV=kVarB_PV(1:nth:end,:);
    kVarC_NoPV=kVarC_NoPV(1:nth:end,:);
    kVarC_PV=kVarC_PV(1:nth:end,:);
    kVar_NoPV=kVar_NoPV(1:nth:end,:);
    kVar_PV=kVar_PV(1:nth:end,:);
    CurA_NoPV=CurA_NoPV(1:nth:end,:);
    CurA_PV=CurA_PV(1:nth:end,:);
    CurB_NoPV=CurB_NoPV(1:nth:end,:);
    CurB_PV=CurB_PV(1:nth:end,:);
    CurC_NoPV=CurC_NoPV(1:nth:end,:);
    CurC_PV=CurC_PV(1:nth:end,:);
else
    nth = 1;
end 

for j=1:15
    SName = SimName;
    %tic
    switch j
        case 1
            Values_DistanceAxis = Dist;
            Values_DataAxis_NoPV = Volt_NoPV;
            Values_DataAxis_WithPV = Volt_PV;
            SName = [SName ' - Voltage'];
        case 2
            Values_DistanceAxis = distA;
            Values_DataAxis_NoPV = VA_NoPV;
            Values_DataAxis_WithPV = VA_PV;
            SName = [SName ' - Voltage (Phase A)'];
        case 3
            Values_DistanceAxis = distB;
            Values_DataAxis_NoPV = VB_NoPV;
            Values_DataAxis_WithPV = VB_PV;
            SName = [SName ' - Voltage (Phase B)'];
        case 4
            Values_DistanceAxis = distC;
            Values_DataAxis_NoPV = VC_NoPV;
            Values_DataAxis_WithPV = VC_PV;
            SName = [SName ' - Voltage (Phase C)'];
        case 5
            Values_DistanceAxis = distA;
            Values_DataAxis_NoPV = kWA_NoPV;
            Values_DataAxis_WithPV = kWA_PV;
            SName = [SName ' - Real Power (Phase A)'];
            Z_Label_3D(1)={'Power, kW'};
            Z_Label_2D(1)={'Power, kW'};
            z_saturation=[0 8000];
            z_range=(-1500:500:1500);
        case 6
            Values_DistanceAxis = distB;
            Values_DataAxis_NoPV = kWB_NoPV;
            Values_DataAxis_WithPV = kWB_PV;
            SName = [SName ' - Real Power (Phase B)'];
        case 7
            Values_DistanceAxis = distC;
            Values_DataAxis_NoPV = kWC_NoPV;
            Values_DataAxis_WithPV = kWC_PV;
            SName = [SName ' - Real Power (Phase C)'];
        case 8
            Values_DistanceAxis = dist;
            Values_DataAxis_NoPV = kW_NoPV;
            Values_DataAxis_WithPV = kW_PV;
            SName = [SName ' - Real Power (Total)'];
            z_range=(-4500:1500:4500);
        case 9
            Values_DistanceAxis = distA;
            Values_DataAxis_NoPV = kVarA_NoPV.*(-1);
            Values_DataAxis_WithPV = kVarA_PV.*(-1);
            SName = [SName ' - Reactive Power (Phase A)'];
            Z_Label_3D(1)={'Power, kVar'};
            Z_Label_2D(1)={'Power, kVar'};
            z_saturation=[-500 4000];
            z_range=(-500:250:1000);
        case 10
            Values_DistanceAxis = distB;
            Values_DataAxis_NoPV = kVarB_NoPV.*(-1);
            Values_DataAxis_WithPV = kVarB_PV.*(-1);
            SName = [SName ' - Reactive Power (Phase B)'];
        case 11
            Values_DistanceAxis = distC;
            Values_DataAxis_NoPV = kVarC_NoPV.*(-1);
            Values_DataAxis_WithPV = kVarC_PV.*(-1);
            SName = [SName ' - Reactive Power (Phase C)'];
        case 12
            Values_DistanceAxis = dist;
            Values_DataAxis_NoPV = kVar_NoPV.*(-1);
            Values_DataAxis_WithPV = kVar_PV.*(-1);
            SName = [SName ' - Reactive Power (Total)'];
            z_range=(-1500:500:1500);
        case 13
            Values_DistanceAxis = distA;
            Values_DataAxis_NoPV = CurA_NoPV;
            Values_DataAxis_WithPV = CurA_PV;
            SName = [SName ' - Current (Phase A)'];
            Z_Label_3D(1)={'Current, Amp'};
            Z_Label_2D(1)={'Current, Amp'};
            z_saturation=[0 550];
            z_range=(-300:150:300);
        case 14
            Values_DistanceAxis = distB;
            Values_DataAxis_NoPV = CurB_NoPV;
            Values_DataAxis_WithPV = CurB_PV;
            SName = [SName ' - Current (Phase B)'];
        case 15
            Values_DistanceAxis = distC;
            Values_DataAxis_NoPV = CurC_NoPV;
            Values_DataAxis_WithPV = CurC_PV;
            SName = [SName ' - Current (Phase C)'];

    end
    
    PlotOpenDSS_3D_Yearly(Values_DistanceAxis,... 
    Values_DataAxis_NoPV,...
    Values_DataAxis_WithPV,...
    Directory_Base_Output,...
    TimeMin,...
    TimeMax,...
    Edge,...
    X_Label_2D,...
    Y_Label_2D,...
    Z_Label_2D,...
    X_Label_3D,...
    Y_Label_3D,...
    Z_Label_3D,...
    Label_FontSize_2D,...
    Label_FontSize_3D,...
    z_saturation,...
    z_range,...
    SName,...
    nth,...
    xrange,...
    DTick)
%%toc
end

%% Output Data
% dssText.Command = 'Export monitor Sub';% Substation monitor
% EventLogEvaluation(dssSolution.EventLog)

end