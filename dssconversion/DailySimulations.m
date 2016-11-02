function [dat] = DailySimulations(Commands, SimName, stepsize, SimPer, p, PVSys, NumPts)

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
% !! Names: this function needs the transformers/capacitors and the regulators/capcontrolers to have
% the same name without any '_' or '.' otherwise EventLogEvaluation
% function won't work!
i=strfind(p,'/');
i=i(length(i));
result_dir=[p(1:i) 'results'];
% result_dir = [p(1:strfind(p,'/f480.dss')) 'results'];
step = str2double(stepsize(1:end-1));

if ~exist(result_dir, 'dir');
    result_dir = [p(1:i) 'results'];
    mkdir(result_dir);
end

dssObj = actxserver('OpendssEngine.dss');
dssText = dssObj.Text;
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
    dssText.Command = ['Loadshape.' PVSys{a} '.Mult = (' num2str(zeros(1,NumPts)+0.00011) ')'];
% checks if the pv loadshapes were set to 0 - just to make sure
%     dssText.Command = ['? Loadshape.' PVSys{a} '.Mult'];
%     dssText.Result
end

dssCircuit = dssObj.ActiveCircuit;
dssSolution = dssCircuit.Solution;

% dssEnergyMeters = dssCircuit.Meters;
dssMonitor = dssCircuit.Monitors;
% Dist = dssCircuit.AllNodeDistances*(3280.8399); % conversion from km to kft (not sure why OpenDSS presents it in km, when the line lengths are clearly identified to be in kft)
% Dist = dssCircuit.AllNodeDistances; % Distance in km
% NodeNames = dssCircuit.AllNodeNames';
dssSolution.MaxControlIterations=1000;
dssSolution.maxiterations=1000;

dssText.Command = 'Set controlmode = static';
dssText.Command = ['Set mode = daily stepsize = ' stepsize];
dssText.Command = 'Set number = 1';% will stop after each solve

dssSolution.InitSnap;
dssSolution.dblHour = 0.0;
 
Volt_MaxMin_NoPV=zeros(NumPts,2); LossTotal_NoPV=zeros(NumPts,2); LossLine_NoPV=zeros(NumPts,2);
Volt_NoPV=zeros(NumPts,length(dssCircuit.AllBusVmagPu)); TotalPower_NoPV=zeros(NumPts,2); 
% Volt_NoPV_PhaseA = []; Volt_NoPV_PhaseB = []; Volt_NoPV_PhaseC = [];
i=1;
nodeName = dssCircuit.AllNodeNames;
%tic
while (dssSolution.dblHour < SimPer)
    dssSolution.Solve;
    if dssSolution.Converged
        Volt_MaxMin_NoPV(i,:) = [max(dssCircuit.AllBusVmagPu(dssCircuit.AllBusVmagPu < 4)) min(dssCircuit.AllBusVmagPu(dssCircuit.AllBusVmagPu > 0.3))];
        LossTotal_NoPV(i,:) = dssCircuit.Losses/(1e+06);
        Volt_NoPV(i,:) = dssCircuit.AllBusVmagPu;
%         Volt_NoPV_PhaseA = [Volt_NoPV_PhaseA; dssCircuit.AllNodeVmagPUByPhase(1)];
%         Volt_NoPV_PhaseB = [Volt_NoPV_PhaseB; dssCircuit.AllNodeVmagPUByPhase(2)];
%         Volt_NoPV_PhaseC = [Volt_NoPV_PhaseC; dssCircuit.AllNodeVmagPUByPhase(3)];
        LossLine_NoPV(i,:) = dssCircuit.LineLosses/(1e+06);
        TotalPower_NoPV(i,:) = dssCircuit.TotalPower/(1e+03)*(-1);
        i=i+1;
    else
        disp(['System did not converge for ''no PV'' case (NICHT GUT) at hour: ' num2str(dssSolution.dblHour)]);
    end
    
end
%toc
circuitElementNames = dssCircuit.AllElementNames;
EventLog = dssSolution.EventLog;
[Event_NoPV Reg_NoPV Cap_NoPV tap2plotNoPV cap2plotNoPV] = EventLogEvaluation(circuitElementNames,EventLog, i-1, step);

disp('- Evaluating EventLog');
CurA_NoPV=zeros(NumPts,2032);  CurB_NoPV=zeros(NumPts,2047); CurC_NoPV=zeros(NumPts,2020);
kWA_NoPV=zeros(NumPts,2032); kWB_NoPV=zeros(NumPts,2047); kWC_NoPV=zeros(NumPts,2020); kW_NoPV=zeros(NumPts,2454);
kVarA_NoPV=zeros(NumPts,2032); kVarB_NoPV=zeros(NumPts,2047); kVarC_NoPV=zeros(NumPts,2020); kVar_NoPV=zeros(NumPts,2454);
% distA=zeros(1,2032); distB=zeros(1,2047); distC=zeros(1,2020); dist=zeros(1,2454);
VA_NoPV=zeros(NumPts,2032); VB_NoPV=zeros(NumPts,2047); VC_NoPV=zeros(NumPts,2020);
i_A=1; i_B=1; i_C=1; i_T=1;
Mon = dssMonitor.First;
%tic

%% for 3D plot
% % while Mon>0,
% %     a = readMonitor(dssMonitor.ByteStream);
% %     data = a.data;DatanoPV=data;
% %     MonName = dssMonitor.Name;
% %     m = regexp(MonName, '\.', 'split');
% %     if length(m)==1
% %         Mon = dssMonitor.next;
% %         continue
% %     end
% %     bus = m{2};
% %     n = cell2mat(m([3:end]));
% %     switch n
% %         case '12'
% %             va=data(:,3);va_ang=data(:,4);vb=data(:,5);vb_ang=data(:,6);
% %             ia=data(:,7);ia_ang=data(:,8);ib=data(:,9);ib_ang=data(:,10);
% %             CurA_NoPV(:,i_A) = ia.*cosd(ia_ang);CurB_NoPV(:,i_B) = ib.*cosd(ib_ang+120);
% %             kWA_NoPV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;kWB_NoPV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;
% %             kW_NoPV(:,i_T) = kWA_NoPV(:,i_A)+kWB_NoPV(:,i_B);
% %             kVarA_NoPV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarB_NoPV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;
% %             kVar_NoPV(:,i_T) = kVarA_NoPV(:,i_A)+kVarB_NoPV(:,i_B);
% %             bus = [bus '.1'];
% %             VA_NoPV(:,i_A) = va./(12000/sqrt(3));VB_NoPV(:,i_B) = vb./(12000/sqrt(3));
% %             i_A=i_A+1; i_B=i_B+1; i_T=i_T+1;
% %         case '13'
% %             va=data(:,3);va_ang=data(:,4);vc=data(:,5);vc_ang=data(:,6);
% %             ia=data(:,7);ia_ang=data(:,8);ic=data(:,9);ic_ang=data(:,10);
% %             CurA_NoPV(:,i_A) = ia.*cosd(ia_ang);CurC_NoPV(:,i_C) = ic.*cosd(ic_ang+120);
% %             kWA_NoPV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;kWC_NoPV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
% %             kW_NoPV(:,i_T) = kWA_NoPV(:,i_A)+kWC_NoPV(:,i_C);
% %             kVarA_NoPV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarC_NoPV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
% %             kVar_NoPV(:,i_T) = kVarA_NoPV(:,i_A)+kVarC_NoPV(:,i_C);
% %             bus = [bus '.1'];
% %             VA_NoPV(:,i_A) = va./(12000/sqrt(3));VC_NoPV(:,i_C) = vc./(12000/sqrt(3));
% %             i_A=i_A+1; i_C=i_C+1; i_T=i_T+1;
% %         case '23'
% %             vb=data(:,3);vb_ang=data(:,4);vc=data(:,5);vc_ang=data(:,6);
% %             ib=data(:,7);ib_ang=data(:,8);ic=data(:,9);ic_ang=data(:,10);
% %             CurB_NoPV(:,i_B) = ib.*cosd(ib_ang);CurC_NoPV(:,i_C) = ic.*cosd(ic_ang+120);
% %             kWB_NoPV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;kWC_NoPV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
% %             kW_NoPV(:,i_T) = kWB_NoPV(:,i_B)+kWC_NoPV(:,i_C);
% %             kVarB_NoPV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;kVarC_NoPV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
% %             kVar_NoPV(:,i_T) = kVarB_NoPV(:,i_B)+kVarC_NoPV(:,i_C);
% %             bus = [bus '.2'];
% %             VB_NoPV(:,i_B) = vb./(12000/sqrt(3));VC_NoPV(:,i_C) = vc./(12000/sqrt(3));
% %             i_B=i_B+1; i_C=i_C+1; i_T=i_T+1;
% %         case '123'
% %             va=data(:,3);va_ang=data(:,4);vb=data(:,5);vb_ang=data(:,6);vc=data(:,7);vc_ang=data(:,8);
% %             ia=data(:,9);ia_ang=data(:,10);ib=data(:,11);ib_ang=data(:,12);ic=data(:,13);ic_ang=data(:,14);
% %             CurA_NoPV(:,i_A) = ia.*cosd(ia_ang);CurB_NoPV(:,i_B) = ib.*cosd(ib_ang+120);CurC_NoPV(:,i_C) = ic.*cosd(ic_ang-120);
% %             kWA_NoPV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;kWB_NoPV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;kWC_NoPV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
% %             kW_NoPV(:,i_T) = kWA_NoPV(:,i_A)+kWB_NoPV(:,i_B)+kWC_NoPV(:,i_C);
% %             kVarA_NoPV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarB_NoPV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;kVarC_NoPV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
% %             kVar_NoPV(:,i_T) = kVarA_NoPV(:,i_A)+kVarB_NoPV(:,i_B)+kVarC_NoPV(:,i_C);
% %             bus = [bus '.1'];
% %             VA_NoPV(:,i_A) = va./(12000/sqrt(3));VB_NoPV(:,i_B) = vb./(12000/sqrt(3));VC_NoPV(:,i_C) = vc./(12000/sqrt(3));
% %             i_A=i_A+1; i_B=i_B+1; i_C=i_C+1; i_T=i_T+1;
% %         case '1'
% %             va=data(:,3);va_ang=data(:,4);
% %             ia=data(:,5);ia_ang=data(:,6);
% %             CurA_NoPV(:,i_A) = ia.*cosd(ia_ang);
% %             kWA_NoPV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;
% %             kW_NoPV(:,i_T) = kWA_NoPV(:,i_A);
% %             kVarA_NoPV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;
% %             kVar_NoPV(:,i_T) = kVarA_NoPV(:,i_A);
% %             bus = [bus '.1'];
% %             VA_NoPV(:,i_A) = va./(12000/sqrt(3));
% %             i_A=i_A+1; i_T=i_T+1;
% %         case '2'
% %             vb=data(:,3);vb_ang=data(:,4);
% %             ib=data(:,5);ib_ang=data(:,6);
% %             CurB_NoPV(:,i_B) = ib.*cosd(ib_ang+120);
% %             kWB_NoPV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;
% %             kW_NoPV(:,i_T) = kWB_NoPV(:,i_B);
% %             kVarB_NoPV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;
% %             kVar_NoPV(:,i_T) = kVarB_NoPV(:,i_B);
% %             bus = [bus '.2'];
% %             VB_NoPV(:,i_B) = vb./(12000/sqrt(3));
% %             i_B=i_B+1; i_T=i_T+1;
% %         case '3'
% %             vc=data(:,3);vc_ang=data(:,4);
% %             ic=data(:,5);ic_ang=data(:,6);
% %             CurC_NoPV(:,i_C) = ic.*cosd(ic_ang-120);
% %             kWC_NoPV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
% %             kW_NoPV(:,i_T) = kWC_NoPV(:,i_C);
% %             kVarC_NoPV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
% %             kVar_NoPV(:,i_T) = kVarC_NoPV(:,i_C);
% %             bus = [bus '.3'];
% %             VC_NoPV(:,i_C) = vc./(12000/sqrt(3));
% %             i_C=i_C+1; i_T=i_T+1;
% %     end
% %     
% %     Mon = dssMonitor.next;
% % end
%toc
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
Volt_PV=zeros(NumPts,length(dssCircuit.AllBusVmagPu)); TotalPower_PV=zeros(NumPts,2);
% Volt_PV_PhaseA = []; Volt_PV_PhaseB = []; Volt_PV_PhaseC = [];
i=1;

%% %% %% %% %% %% TESTING %% %% %% %% %% %% %% %% %% 
% [abc before1] = ismember('04804001.1',nodeName);
% [abc before2] = ismember('04804001.2',nodeName);
% [abc before3] = ismember('04804001.3',nodeName);
% [abc after1] = ismember('048047.1',nodeName);
% [abc after2] = ismember('048047.2',nodeName);
% [abc after3] = ismember('048047.3',nodeName);
% [abc after481] = ismember('048048.1',nodeName);
% [abc after482] = ismember('048048.2',nodeName);
% [abc after483] = ismember('048048.3',nodeName);
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
%tic
for uh = 1:length(dssCircuit.AllElementNames)
    tmpva=dssCircuit.AllElementNames{uh};
    elts(uh) = {tmpva(1:4)};
end
% time of sunset
sunset=19;
%simulation until sunset when it becomes the same as without PV
while (single(dssSolution.dblHour) < SimPer) 
    if ((dssSolution.dblHour==sunset-0.6)|(dssSolution.dblHour==sunset-0.2 ))&~ischar(tap2plotNoPV)
        nmel = find(ismember(elts,'RegC'));
        fieldstap = fields(tap2plotNoPV);
        tapfiname = regexp(fieldstap,'\_','split');
        for kd=1:length(nmel)
            eltsName = regexp(dssCircuit.AllElementNames{nmel(kd)},'\.','split');
            for yh = 1:length(tapfiname)
                tapnames(yh) = lower(tapfiname{yh}(2));
            end
            dssText.Command = [dssCircuit.AllElementNames{nmel(kd)} '.TapNum = ' num2str((tap2plotNoPV.(fieldstap{find(ismember(tapnames,eltsName{2}))})(single(dssSolution.dblHour*2880/24))-1)*32/0.2)];
        end
    end
    if (single(dssSolution.dblHour)<sunset)
        dssSolution.Solve;
        if dssSolution.Converged
            Volt_MaxMin_PV(i,:) = [max(dssCircuit.AllBusVmagPu(dssCircuit.AllBusVmagPu < 4)) min(dssCircuit.AllBusVmagPu(dssCircuit.AllBusVmagPu > 0.3))];
            LossTotal_PV(i,:) = dssCircuit.Losses/(1e+06);
            LossLine_PV(i,:) = dssCircuit.LineLosses/(1e+06);
            Volt_PV(i,:) = dssCircuit.AllBusVmagPu;
    %         Volt_PV_PhaseA = [Volt_PV_PhaseA; dssCircuit.AllNodeVmagPUByPhase(1)];
    %         Volt_PV_PhaseB = [Volt_PV_PhaseB; dssCircuit.AllNodeVmagPUByPhase(2)];
    %         Volt_PV_PhaseC = [Volt_PV_PhaseC; dssCircuit.AllNodeVmagPUByPhase(3)];
            TotalPower_PV(i,:) = dssCircuit.TotalPower/(1e+03)*(-1);
    %         disp(Volt_PV(i,[1070:1072,1397:1399]));
    %         disp(Volt_PV(i,[after1:after3])-Volt_PV(i,[before1:before3]));
    %         disp(Volt_PV(i,[after481:after483])-Volt_PV(i,[before1:before3]));
            i=i+1;
        else
            disp(['System did not converge for ''with PV'' case (NICHT GUT) at hour: ' num2str(dssSolution.dblHour)]);
        end
    else

    Volt_MaxMin_PV(i,:) = Volt_MaxMin_NoPV(i,:);
    LossTotal_PV(i,:) = LossTotal_NoPV(i,:);
    LossLine_PV(i,:) = LossLine_NoPV(i,:);
    Volt_PV(i,:) = Volt_NoPV(i,:);
    TotalPower_PV(i,:) = TotalPower_NoPV(i,:);
    dssSolution.dblHour=dssSolution.dblHour+30/3600;
    i=i+1;
end
    
end
%toc

circuitElementNames = dssCircuit.AllElementNames;
EventLognopv = EventLog;
evt = regexp(EventLognopv,'\,','split');
for kl= 1:length(evt)
    aaa=evt{kl}{1};
    ho=str2num(aaa(6:end));
    aaa=evt{kl}{2};
    se=str2num(aaa(6:end));
    selectEvt(kl) = (sunset<(ho+se/3600));
end
if isempty(evt); selectEvt='';end;
toaddtoEvtPV = EventLognopv(selectEvt);
EventLog = dssSolution.EventLog;
EventLog = vertcat(EventLog,toaddtoEvtPV);
[Event_PV Reg_PV Cap_PV tap2plotPV cap2plotPV] = EventLogEvaluation(circuitElementNames,EventLog, i-1, step); % TO CHANGE
% FOR EVERY CIRCUIT!!! In EventLogEvaluation function all the names of cap
% and reg controlers automatically registered!!


% read all the monitors in the circuit
disp('- Evaluating EventLog');
CurA_PV=zeros(NumPts,2032);  CurB_PV=zeros(NumPts,2047); CurC_PV=zeros(NumPts,2020);
kWA_PV=zeros(NumPts,2032); kWB_PV=zeros(NumPts,2047); kWC_PV=zeros(NumPts,2020); kW_PV=zeros(NumPts,2454);
kVarA_PV=zeros(NumPts,2032); kVarB_PV=zeros(NumPts,2047); kVarC_PV=zeros(NumPts,2020); kVar_PV=zeros(NumPts,2454);
distA=zeros(1,2032); distB=zeros(1,2047); distC=zeros(1,2020); dist=zeros(1,2454);
VA_PV=zeros(NumPts,2032); VB_PV=zeros(NumPts,2047); VC_PV=zeros(NumPts,2020);
i_A=1; i_B=1; i_C=1; i_T=1;
Mon = dssMonitor.First;
%tic

%% for 3D plot
% % while Mon>0,
% %     a = readMonitor(dssMonitor.ByteStream);
% %     data = a.data;data = vertcat(data,DatanoPV(2341:2880,:))
% %     MonName = dssMonitor.Name;
% %     m = regexp(MonName, '\.', 'split');
% %     if length(m)==1
% %         Mon = dssMonitor.next;
% %         continue
% %     end
% %     bus = m{2};
% %     n = cell2mat(m([3:end]));
% %     switch n
% %         case '12'
% %             va=data(:,3);va_ang=data(:,4);vb=data(:,5);vb_ang=data(:,6);
% %             ia=data(:,7);ia_ang=data(:,8);ib=data(:,9);ib_ang=data(:,10);
% %             CurA_PV(:,i_A) = ia.*cosd(ia_ang);CurB_PV(:,i_B) = ib.*cosd(ib_ang+120);
% %             kWA_PV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;kWB_PV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;
% %             kW_PV(:,i_T) = kWA_PV(:,i_A)+kWB_PV(:,i_B);
% %             kVarA_PV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarB_PV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;
% %             kVar_PV(:,i_T) = kVarA_PV(:,i_A)+kVarB_PV(:,i_B);
% %             bus = [bus '.1'];
% %             [TF,LOC] = ismember(bus, NodeNames);
% %             distA(:,i_A) = Dist(LOC); distB(:,i_B) = Dist(LOC); dist(:,i_T) = Dist(LOC);
% %             VA_PV(:,i_A) = va./(12000/sqrt(3)); VB_PV(:,i_B) = vb./(12000/sqrt(3));
% %             i_A=i_A+1; i_B=i_B+1; i_T=i_T+1;
% %         case '13'
% %             va=data(:,3);va_ang=data(:,4);vc=data(:,5);vc_ang=data(:,6);
% %             ia=data(:,7);ia_ang=data(:,8);ic=data(:,9);ic_ang=data(:,10);
% %             CurA_PV(:,i_A) = ia.*cosd(ia_ang);CurC_PV(:,i_C) = ic.*cosd(ic_ang+120);
% %             kWA_PV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;kWC_PV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
% %             kW_PV(:,i_T) = kWA_PV(:,i_A)+kWC_PV(:,i_C);
% %             kVarA_PV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarC_PV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
% %             kVar_PV(:,i_T) = kVarA_PV(:,i_A)+kVarC_PV(:,i_C);
% %             bus = [bus '.1'];
% %             [TF,LOC] = ismember(bus, NodeNames);
% %             distA(:,i_A) = Dist(LOC); distC(:,i_C) = Dist(LOC);dist(:,i_T) = Dist(LOC);
% %             VA_PV(:,i_A) = va./(12000/sqrt(3)); VC_PV(:,i_C) = vc./(12000/sqrt(3));
% %             i_A=i_A+1; i_C=i_C+1; i_T=i_T+1;
% %         case '23'
% %             vb=data(:,3);vb_ang=data(:,4);vc=data(:,5);vc_ang=data(:,6);
% %             ib=data(:,7);ib_ang=data(:,8);ic=data(:,9);ic_ang=data(:,10);
% %             CurB_PV(:,i_B) = ib.*cosd(ib_ang);CurC_PV(:,i_C) = ic.*cosd(ic_ang+120);
% %             kWB_PV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;kWC_PV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
% %             kW_PV(:,i_T) = kWB_PV(:,i_B)+kWC_PV(:,i_C);
% %             kVarB_PV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;kVarC_PV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
% %             kVar_PV(:,i_T) = kVarB_PV(:,i_B)+kVarC_PV(:,i_C);
% %             bus = [bus '.2'];
% %             [TF,LOC] = ismember(bus, NodeNames);
% %             distB(:,i_B) = Dist(LOC); distC(:,i_C) = Dist(LOC); dist(:,i_T) = Dist(LOC);
% %             VB_PV(:,i_B) = vb./(12000/sqrt(3)); VC_PV(:,i_C) = vc./(12000/sqrt(3));
% %             i_B=i_B+1; i_C=i_C+1; i_T=i_T+1;
% %         case '123'
% %             va=data(:,3);va_ang=data(:,4);vb=data(:,5);vb_ang=data(:,6);vc=data(:,7);vc_ang=data(:,8);
% %             ia=data(:,9);ia_ang=data(:,10);ib=data(:,11);ib_ang=data(:,12);ic=data(:,13);ic_ang=data(:,14);
% %             CurA_PV(:,i_A)=ia.*cosd(ia_ang); CurB_PV(:,i_B) = ib.*cosd(ib_ang+120);CurC_PV(:,i_C) = ic.*cosd(ic_ang-120);
% %             kWA_PV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)./1000;kWB_PV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;kWC_PV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
% %             kW_PV(:,i_T) = kWA_PV(:,i_A)+kWB_PV(:,i_B)+kWC_PV(:,i_C);
% %             kVarA_PV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;kVarB_PV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;kVarC_PV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
% %             kVar_PV(:,i_T) = kVarA_PV(:,i_A)+kVarB_PV(:,i_B)+kVarC_PV(:,i_C);
% %             bus = [bus '.1'];
% %             [TF,LOC] = ismember(bus, NodeNames);
% %             distA(:,i_A) = Dist(LOC); distB(:,i_B) = Dist(LOC); distC(:,i_C) = Dist(LOC);dist(:,i_T) = Dist(LOC);
% %             VA_PV(:,i_A) = va./(12000/sqrt(3)); VB_PV(:,i_B) = vb./(12000/sqrt(3)); VC_PV(:,i_C) = vc./(12000/sqrt(3));
% %             i_A=i_A+1; i_B=i_B+1; i_C=i_C+1; i_T=i_T+1;
% %         case '1'
% %             va=data(:,3);va_ang=data(:,4);
% %             ia=data(:,5);ia_ang=data(:,6);
% %             CurA_PV(:,i_A) = ia.*cosd(ia_ang);
% %             kWA_PV(:,i_A) = va.*ia.*cosd(va_ang-ia_ang)/1000;
% %             kW_PV(:,i_T) = kWA_PV(:,i_A);
% %             kVarA_PV(:,i_A) = va.*ia.*sind(va_ang-ia_ang)/1000;
% %             kVar_PV(:,i_T) = kVarA_PV(:,i_A);
% %             bus = [bus '.1'];
% %             [TF,LOC] = ismember(bus, NodeNames);
% %             distA(:,i_A) = Dist(LOC); dist(:,i_T) = Dist(LOC);
% %             VA_PV(:,i_A) = va./(12000/sqrt(3));
% %             i_A=i_A+1;  i_T=i_T+1;
% %         case '2'
% %             vb=data(:,3);vb_ang=data(:,4);
% %             ib=data(:,5);ib_ang=data(:,6);
% %             CurB_PV(:,i_B) = ib.*cosd(ib_ang+120);
% %             kWB_PV(:,i_B) = vb.*ib.*cosd(vb_ang-ib_ang)/1000;
% %             kW_PV(:,i_T) = kWB_PV(:,i_B);
% %             kVarB_PV(:,i_B) = vb.*ib.*sind(vb_ang-ib_ang)/1000;
% %             kVar_PV(:,i_T) = kVarB_PV(:,i_B);
% %             bus = [bus '.2'];
% %             [TF,LOC] = ismember(bus, NodeNames);
% %             distB(:,i_B) = Dist(LOC); dist(:,i_T) = Dist(LOC);
% %             VB_PV(:,i_B) = vb./(12000/sqrt(3));
% %             i_B=i_B+1;  i_T=i_T+1;
% %         case '3'
% %             vc=data(:,3);vc_ang=data(:,4);
% %             ic=data(:,5);ic_ang=data(:,6);
% %             CurC_PV(:,i_C) = ic.*cosd(ic_ang-120);
% %             kWC_PV(:,i_C) = vc.*ic.*cosd(vc_ang-ic_ang)/1000;
% %             kW_PV(:,i_T) = kWC_PV(:,i_C);
% %             kVarC_PV(:,i_C) = vc.*ic.*sind(vc_ang-ic_ang)/1000;
% %             kVar_PV(:,i_T) = kVarC_PV(:,i_C);
% %             bus = [bus '.3'];
% %             [TF,LOC] = ismember(bus, NodeNames);
% %             distC(:,i_C) = Dist(LOC);dist(:,i_T) = Dist(LOC);
% %             VC_PV(:,i_C) = vc./(12000/sqrt(3));
% %             i_C=i_C+1;  i_T=i_T+1;
% %     end
% %     
% %     Mon = dssMonitor.next;
% % end
%toc
% OPTIONAL: exports all monitors to hard drive. See save monitors section
% in main_f520.m
% dssText.Command = 'cd C:\Work\Projects\2012\1787-UCSD_PV\Simulation\System\SDGE\tmp\f520\Monitors';
% dssText.Command = 'redirect f520_mon.dss';

%%
if ~isempty(Reg_PV)
    RegTotal_PV = sum(Reg_PV(2:end,end));
    RegTotal_NoPV = sum(Reg_NoPV(2:end,end));
end
if ~isempty(Cap_PV)
    CapTotal_PV = sum(Cap_PV(2:end,end));
    CapTotal_NoPV = sum(Cap_NoPV(2:end,end));
end
if exist('RegTotal_PV'); 
    Operations.RegTotal_PV = RegTotal_PV;
    Operations.RegTotal_NoPV = RegTotal_NoPV;
else 
    Operations.RegTotal_PV = 0;
    Operations.RegTotal_NoPV = 0;
end
if exist('CapTotal_PV');
    Operations.CapTotal_PV = CapTotal_PV;
    Operations.CapTotal_NoPV = CapTotal_NoPV;
else
    Operations.CapTotal_PV = 0;
    Operations.CapTotal_NoPV = 0;
end
save([result_dir '/' SimName '_results.mat'], 'LossTotal_PV', 'LossTotal_NoPV', 'TotalPower_PV', 'TotalPower_NoPV','Volt_PV','nodeName','Volt_NoPV','Operations');

%% 2D Plots
time = i-1;
xrange = 1:(i-1);
if NumPts > 2000
    time=time/120;
    xrange = xrange./120;
end
colorCurve = [{'r*-'} , {'g*-'},{'b*-'},{'k*-'},{'y*-'},{'m*-'},{'rx:'},{'gx:'},{'bx:'},{'kx:'}];
ColorLine = [{'k'} , {'b'},{'c'},{'y'},{'m'},{'r'},{'g'}];
disp(['Plotting results for ' SimName]);
%tic
% if the visibility for the 2D plots is set to 'off', Matlab does a poor
% job saving .fig files - they don't open for editing, anymore. the
% simulaiton gets slower but at least more editing can be done later on.
% f = figure('visible','off'); hold on;
f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
h = plot(xrange, Volt_MaxMin_PV(:,1), 'r', xrange, Volt_MaxMin_NoPV(:,1),...
    'g', xrange, Volt_MaxMin_PV(:,2), 'b', xrange, Volt_MaxMin_NoPV(:,2), 'k');
legend(h,{'Max (PV)', 'Max (No PV)','Min (PV)', 'Min (No PV)'},'fontsize',20);
title('Maximum and Minimum Bus Voltages','fontsize',20)
xlabel('Time of Day','fontsize',20)
ylabel('Voltages (pu)','fontsize',20)
grid on;
xlim([min(xrange) time])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')


% f = figure('units','normalized','outerposition',[0 0 1 1])('visible','off'); hold on;
f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
h = plot(xrange, LossTotal_PV(:,1), 'r', xrange, LossTotal_NoPV(:,1),...
    'g', xrange, LossTotal_PV(:,2), 'b', xrange, LossTotal_NoPV(:,2), 'k');
legend(h,{'MW (PV)', 'MW (No PV)', 'MVar (PV)', 'MVar (No PV)'},'fontsize',20)
title('Total Losses','fontsize',20)
xlabel('Time of Day','fontsize',20)
xlim([min(xrange) time])
ylabel('Losses','fontsize',20)
grid on;
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

% f = figure('visible','off'); hold on;
f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
h = plot(xrange, LossLine_PV(:,1), 'r', xrange, LossLine_NoPV(:,1),...
    'g', xrange, LossLine_PV(:,2), 'b', xrange, LossLine_NoPV(:,2), 'k');
legend(h,{'MW (PV)', 'MW (No PV)', 'MVar (PV)', 'MVar (No PV)'},'fontsize',20)
title('Total Line Losses','fontsize',20)
xlabel('Time of Day','fontsize',20)
xlim([min(xrange) time])
ylabel('Losses','fontsize',20)
grid on;
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

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
if ~isempty(Event_PV.RegControl)
    for ij = 1:length(Names_Reg)
        RegCtr = Names_Reg{ij};
        for ii=2:length(AA.(RegCtr))
            line(ones(1,2)*AA.(RegCtr)(ii)/3600,[min(min(min([TotalPower_PV(:,1)]),min([TotalPower_PV(:,2)]))) max(max([TotalPower_PV(:,1)]))],'Color',ColorLine{ij},'LineStyle','--');
            if length(AA.(RegCtr))>2
                if ii==3
                text(AA.(RegCtr)(ii)/3600,min(min(min([TotalPower_PV(:,1)]),min([TotalPower_PV(:,2)]))),['\leftarrow Tap op. ' RegCtr],'FontSize',16)
                end
            end
        end
    end
end
    legend(h,{'MW (PV)', 'MW (No PV)', 'MVar (PV)', 'MVar (No PV)'},'fontsize',20)
    title('Total Power','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('Power','fontsize',20)
    grid on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

close all;
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
            grid on;
        end
    catch
        for z=1:numOfRegCurve
            plot(xrange, Reg_PV((2:length(Reg_PV)/length(xrange):end),z),colorCurve{z});
            grid on;
        end
    end
    legend([LegendList{:}])
    title('Voltage Regulator Events (With PV)','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('# Events','fontsize',20)
    grid on;
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
            grid on;
        end
    catch
        for z=1:numOfRegCurve
            plot(xrange, Reg_NoPV((2:length(Reg_NoPV)/length(xrange):end),z),colorCurve{z});
            grid on;
        end
    end
    legend([LegendList{:}])
    title('Voltage Regulator Events (No PV)','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('# Events','fontsize',20)
    grid on;
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
        grid on;
    end
    legend([LegendList{:}])
    title('Capacitor Bank Events (With PV)','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('# Events','fontsize',20)
    grid on;
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
        grid on;
    end
    legend([LegendList{:}])
    title('Capacitor Bank Events (No PV)','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('# Events','fontsize',20)
    grid on;
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
    legend(h,{['Total_P_V (' num2str(RegTotal_PV) ')'], ['Total_N_o_P_V (' num2str(RegTotal_NoPV) ')']})
    title('Voltage Regulator Events - Total','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('# Events','fontsize',20)
    grid on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
end
if ~isempty(Cap_PV)
    CapTotal_PV = sum(Cap_PV(2:end,end));
    CapTotal_NoPV = sum(Cap_NoPV(2:end,end));
    % f = figure('visible','off'); hold on;
    f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    h = plot(xrange, Cap_PV(2:end,end), 'r', xrange, Cap_NoPV(2:end,end), 'g');
    legend(h,{['Total_P_V (' num2str(CapTotal_PV) ')'], ['Total_N_o_P_V (' num2str(CapTotal_NoPV) ')']})
    title('Capacitor Bank Events - Total','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('# Events','fontsize',20)
    grid on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
end


%% Plot specific node to Point Loma circuit
% Voltage at substation and end of lines - Phase 1
%auto code:
subBus = dssCircuit.AllBusNames(2);
endBus = dssCircuit.AllBusNames(find(ismember(dssCircuit.AllBusDistances,max(dssCircuit.AllBusDistances))));
[abc sub1] = ismember([subBus{:} '.1'],nodeName);
[abc end_1] = ismember([endBus{:} '.1'],nodeName);
f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
h = plot(xrange, Volt_PV(:,sub1), 'r', xrange, Volt_PV(:,end_1), 'g', xrange, Volt_PV(:,end_1), 'b');
legend(h,{['Voltage at substation - P1'], ['Voltage at end of line ' endBus{:} ' - P1']})
title('Voltage at substation and the furthest end of line - Phase 1','fontsize',20)
xlabel('Time of Day','fontsize',20)
xlim([min(xrange) time])
ylabel('Voltage (pu)','fontsize',20)
grid on;
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
  
% % [abc sub1] = ismember('0480.1',nodeName);
% % [abc end4519_1] = ismember('04804519.1',nodeName);
% % [abc end4751_1] = ismember('04804751.1',nodeName);
% % [abc end9926_1] = ismember('04809926.3',nodeName);
% % 
% % f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
% % h = plot(xrange, Volt_PV(:,sub1), 'r', xrange, Volt_PV(:,end4519_1), 'g', xrange, Volt_PV(:,end4751_1), 'b');
% % legend(h,{['Voltage at substation - P1'], ['Voltage at end of line 04804519 - P1'],['Voltage at end of line 04804751 - P1']})
% % title('Voltage at substation and 2 main ends of line - Phase 1','fontsize',20)
% % xlabel('Time of Day','fontsize',20)
% % xlim([min(xrange) time])
% % ylabel('Voltage (pu)','fontsize',20)
% % grid on;
% % saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
% % saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

% Voltage at substation and end of lines - Phase 2
if ismember('048043.1',nodeName)
    [abc sub2] = ismember('0480.2',nodeName);
    [abc end4519_2] = ismember('04804519.2',nodeName);
    [abc end4751_2] = ismember('04804751.2',nodeName);
    [abc end9926_2] = ismember('04809926.2',nodeName);

    f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    h = plot(xrange, Volt_PV(:,sub2), 'r', xrange, Volt_PV(:,end4519_2), 'g', xrange, Volt_PV(:,end4751_2), 'b',xrange, Volt_PV(:,end9926_2), 'm');
    legend(h,{['Voltage at substation - P2'], ['Voltage at end of line 04804519 - P2'],['Voltage at end of line 04804751 - P2'],['Voltage at end of line 04809926 - P2']})
    title('Voltage at substation and 3 main ends of line - Phase 2','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('Voltage (pu)','fontsize',20)
    grid on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

    % Voltage at substation and end of lines - Phase 3

    [abc sub3] = ismember('0480.3',nodeName);
    [abc end4519_3] = ismember('04804519.3',nodeName);
    [abc end4751_3] = ismember('04804751.3',nodeName);
    [abc end9926_3] = ismember('04809926.3',nodeName);

    f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    h = plot(xrange, Volt_PV(:,sub3), 'r', xrange, Volt_PV(:,end4519_3), 'g', xrange, Volt_PV(:,end4751_3), 'b',xrange, Volt_PV(:,end9926_3), 'm');
    legend(h,{['Voltage at substation - P3'], ['Voltage at end of line 04804519 - P3'],['Voltage at end of line 04804751 - P3'],['Voltage at end of line 04809926 - P3']})
    title('Voltage at substation and 3 main ends of line - Phase 3','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('Voltage (pu)','fontsize',20)
    grid on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
    close all;
end

%Code to make it automatic for other circuits problem: dssparse is quite
%long because the circuit has a lot of devices.
close all;
disp('dssparse running');tic;
MatlabCircuit = dssparse(p);
toc;
if isfield(MatlabCircuit,'regcontrol')
    for numRegControlers = 1:length(MatlabCircuit.regcontrol)
        numTransf = find(ismember(lower({MatlabCircuit.transformer.Name}),lower(MatlabCircuit.regcontrol(numRegControlers).transformer)));
        Buses = MatlabCircuit.transformer(numTransf).Buses;
        Buses(1) = regexp(Buses(1),'\.','split'); 
        asds= Buses(1);
        Buses(1)=asds{1}(1);
        Buses(2) = regexp(Buses(2),'\.','split'); 
        asds= Buses(2);
        Buses(2)=asds{1}(1);
        [abc before1] = ismember([lower(Buses{1}) '.1'],nodeName);
        [abc before2] = ismember([lower(Buses{1}) '.2'],nodeName);
        [abc before3] = ismember([lower(Buses{1}) '.3'],nodeName);
        [abc after1] = ismember([lower(Buses{2}) '.1'],nodeName);
        [abc after2] = ismember([lower(Buses{2}) '.2'],nodeName);
        [abc after3] = ismember([lower(Buses{2}) '.3'],nodeName);
        if ~isempty(Event_PV.RegControl) & isfield(Event_PV.RegControl,['reg_' lower(MatlabCircuit.regcontrol(numRegControlers).transformer)])
            regName = ['reg_' lower(MatlabCircuit.regcontrol(numRegControlers).transformer)];
            AA=unique([Event_PV.RegControl.(regName).Regulations.time]);
            f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
            h = plot(xrange, Volt_PV(:,before1), 'r', ...%xrange, Volt_PV(:,before2), 'm', xrange, Volt_PV(:,before3), 'g',...
            xrange, Volt_PV(:,after1), 'k');%, xrange, Volt_PV(:,after2), 'b', xrange, Volt_PV(:,after3), 'c');
            for ii=2:length(AA)
            line(ones(1,2)*AA(ii)/3600,[min(min([Volt_PV(:,before1)]-0.1)) max(max([Volt_PV(:,before1)]-0.1))],'Color','r','LineStyle','--');
            end
            legend(h,{['Voltage before VReg 1 (Phase 1)'], ...['Voltage before VReg 1 (Phase 2)'],['Voltage before VReg 1 (Phase 3)'],...
            ['Voltage after VReg 1 (Phase 1)']});% , ['Voltage after VReg 1 (Phase 2)'],['Voltage after VReg 1 (Phase 3)']})
            title('Voltage before and after Vreg 1','fontsize',20)
            xlabel('Time of Day','fontsize',20)
            xlim([min(xrange) time])
            ylabel('Voltage (pu)','fontsize',20)
            grid on;
            saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
            saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
        end
    end
end
%Biggest load
biggestloadBus = MatlabCircuit.load(find(ismember([MatlabCircuit.load.kVA],max([MatlabCircuit.load.kVA])))).bus1;
if iscell(biggestloadBus);  biggestloadBus = biggestloadBus{1};end
biggestloadBus = regexp(biggestloadBus,'\.','split');
biggestloadBus=biggestloadBus{1};
biggestload1 = find(ismember(nodeName,[biggestloadBus '.1']));
biggestload2 = find(ismember(nodeName,[biggestloadBus '.2']));
biggestload3 = find(ismember(nodeName,[biggestloadBus '.3']));
biggestload = {biggestload1 biggestload2 biggestload3};
colorsplot = {'r' 'k' 'm' 'b' 'g' 'c'};
leg = {['Voltage at biggest load with PV (Phase 1)'] , ['Voltage at biggest load without PV (Phase 1)'] , ...
    ['Voltage at biggest load with PV (Phase 2)'] , ['Voltage at biggest load without PV (Phase 2)'] , ...
    ['Voltage at biggest load with PV (Phase 3)'] , ['Voltage at biggest load without PV (Phase 3)']};
f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for pl=1:3
    try
        plot(xrange, Volt_PV(:,biggestload{pl}),colorsplot{2*pl-1});
        plot(xrange, Volt_NoPV(:,biggestload{pl}),colorsplot{2*pl});
    catch err
        leg(2*pl)='';
        leg(2*pl-1)='';
    end
end
legend(leg)
title('Voltage at biggest load','fontsize',20)
xlabel('Time of Day','fontsize',20)
xlim([min(xrange) time])
ylabel('Voltage (pu)','fontsize',20)
grid on;
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
%
%auto    
% %     [abc before1] = ismember('048043.1',nodeName);
% %     [abc before2] = ismember('048043.2',nodeName);
% %     [abc before3] = ismember('048043.3',nodeName);
% %     [abc after1] = ismember('048044.1',nodeName);
% %     [abc after2] = ismember('048044.2',nodeName);
% %     [abc after3] = ismember('048044.3',nodeName);
% %     if ~isempty(Event_PV.RegControl) & isfield(Event_PV.RegControl,'reg_t3')
% %         AA=unique([Event_PV.RegControl.reg_t3.Regulations.time]);
% %         f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
% %         h = plot(xrange, Volt_PV(:,before1), 'r', ...%xrange, Volt_PV(:,before2), 'm', xrange, Volt_PV(:,before3), 'g',...
% %             xrange, Volt_PV(:,after1), 'k');%, xrange, Volt_PV(:,after2), 'b', xrange, Volt_PV(:,after3), 'c');
% %         for ii=2:length(AA)
% %             line(ones(1,2)*AA(ii)/3600,[min(min([Volt_PV(:,before1)]-0.1)) max(max([Volt_PV(:,before1)]-0.1))],'Color','r','LineStyle','--');
% %         end
% %         legend(h,{['Voltage before VReg 1 (Phase 1)'], ...['Voltage before VReg 1 (Phase 2)'],['Voltage before VReg 1 (Phase 3)'],...
% %             ['Voltage after VReg 1 (Phase 1)']});% , ['Voltage after VReg 1 (Phase 2)'],['Voltage after VReg 1 (Phase 3)']})
% %         title('Voltage before and after Vreg 1','fontsize',20)
% %         xlabel('Time of Day','fontsize',20)
% %         xlim([min(xrange) time])
% %         ylabel('Voltage (pu)','fontsize',20)
% %         grid on;
% %         saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
% %         saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
% %     end
% % %auto
% %     if ~isempty(Event_PV.RegControl) & isfield(Event_PV.RegControl,'reg_t4')
% %         [abc before1] = ismember('048056.1',nodeName);
% %         [abc before2] = ismember('048056.2',nodeName);
% %         [abc before3] = ismember('048056.3',nodeName);
% %         [abc after1] = ismember('048057.1',nodeName);
% %         [abc after2] = ismember('048057.2',nodeName);
% %         [abc after3] = ismember('048057.3',nodeName);
% %         AA=unique([Event_PV.RegControl.reg_t4.Regulations.time]);
% %         f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
% %         h = plot(xrange, Volt_PV(:,before1), 'r',... % xrange, Volt_PV(:,before2), 'm', xrange, Volt_PV(:,before3), 'g',...
% %             xrange, Volt_PV(:,after1), 'k'); %, xrange, Volt_PV(:,after2), 'b', xrange, Volt_PV(:,after3), 'c');
% %         for ii=2:length(AA)
% %             line(ones(1,2)*AA(ii)/3600,[min(min([Volt_PV(:,before1)])) max(max([Volt_PV(:,before1)]))],'Color','r','LineStyle','--');
% %         end
% %         legend(h,{['Voltage before VReg 2 (Phase 1)'],... ['Voltage before VReg 2 (Phase 2)'],['Voltage before VReg 2 (Phase 3)'],...
% %             ['Voltage after VReg 2 (Phase 1)']}); %, ['Voltage after VReg 2 (Phase 2)'],['Voltage after VReg 2 (Phase 3)']})
% %         title('Voltage before and after Vreg 2','fontsize',20)
% %         xlabel('Time of Day','fontsize',20)
% %         xlim([min(xrange) time])
% %         ylabel('Voltage (pu)','fontsize',20)
% %         grid on;
% %         saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
% %         saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
% %     end
if ismember('048043.1',nodeName)
    %78/522 5/215.9 29/2.94
    smallone1 = find(ismember(nodeName,'04809108.1'));
    smallone2 = find(ismember(nodeName,'04809108.2'));
    mediumone1 = find(ismember(nodeName,'04804536.1'));
    mediumone2 = find(ismember(nodeName,'04804536.2'));
    mediumone3 = find(ismember(nodeName,'04804536.3'));
    biggestone1 = find(ismember(nodeName,'04804506.1'));
    biggestone2 = find(ismember(nodeName,'04804506.2'));
    biggestone3 = find(ismember(nodeName,'04804506.3'));
    f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    h = plot(xrange, Volt_PV(:,smallone1), 'r', xrange, Volt_PV(:,smallone2), 'm',...
        xrange, Volt_NoPV(:,smallone1), 'k', xrange, Volt_NoPV(:,smallone2), 'b');
    legend(h,{['Voltage Small PV system (Phase 1)'], ['Voltage Small PV system (Phase 2)'],...
        ['Voltage without PV (Phase 1)'], ['Voltage without PV (Phase 2)']})
    title('Voltage with and without PV at a small PV system','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('Voltage (pu)','fontsize',20)
    grid on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

    f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    h = plot(xrange, Volt_PV(:,mediumone1), 'r', xrange, Volt_PV(:,mediumone3), 'm',...
        xrange, Volt_NoPV(:,mediumone1), 'k', xrange, Volt_NoPV(:,mediumone3), 'b');
    legend(h,{['Voltage Medium PV system (Phase 1)'], ['Voltage Medium PV system (Phase 2)'],...
        ['Voltage without PV (Phase 1)'], ['Voltage without PV (Phase 2)']})
    title('Voltage with and without PV at a Medium PV system','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('Voltage (pu)','fontsize',20)
    grid on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')

    f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    h = plot(xrange, Volt_PV(:,biggestone1), 'r', xrange, Volt_PV(:,biggestone2), 'm', xrange, Volt_PV(:,biggestone3), 'g',...
        xrange, Volt_NoPV(:,biggestone1), 'k', xrange, Volt_NoPV(:,biggestone2), 'b', xrange, Volt_NoPV(:,biggestone3), 'c');
    legend(h,{['Voltage Big PV system (Phase 1)'], ['Voltage Big PV system (Phase 2)'],['Voltage Big PV system (Phase 3)'],...
        ['Voltage without PV (Phase 1)'], ['Voltage without PV (Phase 2)'],['Voltage without PV (Phase 3)']})
    title('Voltage with and without PV at a Big PV system','fontsize',20)
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('Voltage (pu)','fontsize',20)
    grid on;
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
%auto
% %     biggestload1 = find(ismember(nodeName,'04804487.1'));
% %     biggestload2 = find(ismember(nodeName,'04804487.2'));
% %     biggestload3 = find(ismember(nodeName,'04804487.3'));
% %     
% %     f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
% %     h = plot(xrange, Volt_PV(:,biggestload1), 'r', xrange, Volt_PV(:,biggestload2), 'm', xrange, Volt_PV(:,biggestload3), 'g',...
% %         xrange, Volt_NoPV(:,biggestload1), 'k', xrange, Volt_NoPV(:,biggestload2), 'b', xrange, Volt_NoPV(:,biggestload3), 'c');
% %     legend(h,{['Voltage at biggest load with PV (Phase 1)'], ['Voltage at biggest load with PV (Phase 2)'],['Voltage at biggest load with PV (Phase 3)'],...
% %         ['Voltage at biggest load without PV (Phase 1)'], ['Voltage at biggest load without PV (Phase 2)'],['Voltage at biggest load without PV (Phase 3)']})
% %     title('Voltage at biggest load','fontsize',20)
% %     xlabel('Time of Day','fontsize',20)
% %     xlim([min(xrange) time])
% %     ylabel('Voltage (pu)','fontsize',20)
% %     grid on;
% %     saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
% %     saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
end
close all;
%% Tap operation plots
if ~isempty(tap2plotPV)
    RegName = fieldnames(tap2plotPV);
    for nbReg=1:length(RegName)
        f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        h = plot(xrange,[tap2plotPV.(RegName{nbReg})],'b',xrange,[tap2plotNoPV.(RegName{nbReg})],'r');
        legend(h,{['With PV'],['Without PV']})
        title(['Voltage Regulator ' num2str(nbReg) ' (' RegName{nbReg} ') Tap positions'],'fontsize',20)
        xlabel('Time of Day','fontsize',20)
        xlim([min(xrange) time])
        ylabel('(pu)','fontsize',20)
        grid on;
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
        grid on;
    end
    legend([LegendList{:}],'fontsize',20)
    grid on;
    xlabel('Time of Day','fontsize',20)
    xlim([min(xrange) time])
    ylabel('Voltage Regulators Tap positions (pu)','fontsize',20)
    subplot(2,1,1,'align'), plot(xrange, TotalPower_PV(:,1), 'r', xrange, TotalPower_NoPV(:,1),...
        'g', xrange, TotalPower_PV(:,2), 'b', xrange, TotalPower_NoPV(:,2), 'k');
    legend({'MW (PV)', 'MW (No PV)', 'MVar (PV)', 'MVar (No PV)'},'fontsize',20)
    title('Total Power compare to tap position','fontsize',20)
    xlim([min(xrange) time])
    grid on;
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
        legend(h,{['With PV'],['Without PV']})
        title(['Capacitor Bank ' num2str(nbCap) ' (' CapName{nbCap} ') step Positions'],'fontsize',20)
        xlabel('Time of Day','fontsize',20)
        xlim([min(xrange) time])
        ylabel('Step','fontsize',20)
        grid on;
        saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
        saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
    end
end

close all;
%% Plot voltage along the circuit
% change it to plot or not voltage along the feeder at each hour. *****
if 1==1
    volt_pv_plot = Volt_PV;
    volt_pv_plot([2<Volt_NoPV]) = nan;
    volt_pv_plot([0.5>Volt_NoPV]) = nan;
    volt_nopv_plot = Volt_NoPV;
    volt_nopv_plot([2<Volt_NoPV]) = nan;
    volt_nopv_plot([0.5>Volt_NoPV]) = nan;
    %% pb voltage is ..x.. 
    dist_plot = Dist(~isnan(volt_pv_plot(1,:)));
    [order order] = sort(dist_plot);
    %% Plot voltage 
    for time = 120:120:2880;
        f = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        for i=1:length(dist_plot)
            plot([dist_plot(order(i)) dist_plot(order(i))],[volt_pv_plot(time,order(i)) volt_nopv_plot(time,order(i))] ,'r');
        end
        h = plot(dist_plot(order),volt_pv_plot(time,order),'x',dist_plot(order),volt_nopv_plot(time,order),'.');
        legend(h,{'With PV','Without PV'},'fontsize',20);
        grid on;
        xlim([0 max(dist_plot)])
        xlabel('Distance, km','fontsize',20);
        ylabel('Voltage, pu','fontsize',20);
        title(['Voltage at ' num2str(time/120)],'fontsize',20)
        saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.png'])
    %     saveas(f,[result_dir '/' SimName '_' get(get(gca,'title'),'string'), '.fig'], 'fig')
        if (time==120*8) || (time==120*16)
            close all;
        end
    end
    
end
close all
sumpv = zeros(1,2880)
for k=1:length(MatlabCircuit.pvsystem)
    pro = find(ismember({MatlabCircuit.loadshape.Name},MatlabCircuit.pvsystem(k).daily));
    sumpv=MatlabCircuit.loadshape(pro).Mult.*MatlabCircuit.pvsystem(k).kVA+sumpv;
end
figure,
plot(xrange,TotalPower_NoPV(:,1),'r',xrange,TotalPower_PV(:,1)'+LossTotal_NoPV(:,1)'-LossTotal_PV(:,1)'+sumpv/1000,'b')

%toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D Plots
if 1==0
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

    z_saturation=[0.90 1.15];% z_saturation influences color spectrum, the higher the range, the farther away from extreme color
    z_range=(0.9:0.05:1.10); % range of values displayed, also determines z-axis numbering

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
                z_saturation=[0 11000];
                z_range=(-500:500:3000);
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
                z_range=(-1000:3000:8000);
            case 9
                Values_DistanceAxis = distA;
                Values_DataAxis_NoPV = kVarA_NoPV.*(-1);
                Values_DataAxis_WithPV = kVarA_PV.*(-1);
                SName = [SName ' - Reactive Power (Phase A)'];
                Z_Label_3D(1)={'Power, kVar'};
                Z_Label_2D(1)={'Power, kVar'};
                z_saturation=[-500 4000];
                z_range=(-500:500:1500);
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
                z_range=(-500:900:4000);
            case 13
                Values_DistanceAxis = distA;
                Values_DataAxis_NoPV = CurA_NoPV;
                Values_DataAxis_WithPV = CurA_PV;
                SName = [SName ' - Current (Phase A)'];
                Z_Label_3D(1)={'Current, Amp'};
                Z_Label_2D(1)={'Current, Amp'};
                z_saturation=[0 550];
                z_range=(-50:50:400);
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

        PlotOpenDSS_3D(Values_DistanceAxis,... 
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
        nth)
    %%toc
    end
end
%% Output Data
% dssText.Command = 'Export monitor Sub';% Substation monitor
% EventLogEvaluation(dssSolution.EventLog)

end