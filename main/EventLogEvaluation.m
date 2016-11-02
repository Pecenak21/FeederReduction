function [EventLog, Reg, Cap, tap2plot, cap2plot] = EventLogEvaluation(circuit,Event, Pts, step)

%evaluated the EventLog for the current OpenDSS simulation.
%INPUTS:
%"Event" - event log (capacitor switching, regulator tap changes)
%"Pts" - number of 'solves' (e.g. 24 for a daily
%simulaiton with 1h stepsize)
%"circuit" - circuit to create a cap and reg control list

%OUTPUTS:
%"EventLog" - structure containing all the events sorted by control devices,
%hour of simulaiton, and performed action
%"Reg" - array summarizing RegControl events. first row represents the total
%events for the corresponding controller. The remaining rows represent the
%point of time in the current simulation. The columns represent the
%individual controllers lined up in the following manner:
% 'I' 'II' 'III' 'IV' 'V' 'VI' 'Sub'.
%"Cap" - array summarizing CapControl events (same structure as "Reg").
%"tap2plot" - tap position of Voltage regulators
%"tap2plot" - cap position of capcontrolers

%disp('\n\nCreating EventLog...');
EventLog = struct('RegControl', {''}, 'CapControl', {''});

%% Create RegControl and CapControl lists
if isfield(circuit,'regcontrol')
    for k=1:length(circuit.transformer)
        EventLog.RegControl.(fnSanitize(lower(circuit.transformer(k).Name),'x')) = struct('ID', circuit.transformer(k).Name, 'Regulations', {''});
        EventLog.RegControl.(fnSanitize(lower(circuit.transformer(k).Name),'x')).Regulations = struct('time', [], 'action', '');
    end
else
    Reg='';
end
if isfield(circuit,'capcontrol')
    for k=1:length(circuit.capacitor)
        EventLog.CapControl.(fnSanitize(lower(circuit.capacitor(k).Name),'x')) = struct('ID', circuit.capacitor(k).Name, 'Regulations', {''});
        EventLog.CapControl.(fnSanitize(lower(circuit.capacitor(k).Name),'x')).Regulations = struct('time', [], 'action', '');;
    end
else
    Cap='';
    cap2plot = '';
end
%% Sort Events
for i=1:length(Event)
    AA = textscan(Event{i},'%s','delimiter',','); AA = AA{1};
    if strcmpi(AA{4}(9),'R')
        Name = fnSanitize(lower(AA{4}(19:end)),'x');
        if ~isfield(EventLog.RegControl, Name)
            error(['The RegControl ' Name ' is not on the current list. Check if the transformer and the regcontrol have the same name.']);
        end
        Hour = AA{1}(6:end);
        Sec = AA{2}(5:end);
        Action = AA{5}(9:end-1);
        if isempty(EventLog.RegControl.(Name).Regulations(1).time)
            rid = 1;
        else
            rid = length(EventLog.RegControl.(Name).Regulations) + 1;
        end
        EventLog.RegControl.(Name).Regulations(rid).time = (str2double(Hour)*3600)+str2double(Sec);
        EventLog.RegControl.(Name).Regulations(rid).action = Action;
    else
        Name = fnSanitize(lower(AA{4}(19:end)),'x');
        %         Name=['cap_' Name];
        if ~isfield(EventLog.CapControl, Name)
            error(['The CapControl ' Name ' is not on the current list. Please add it in!']);
        end
        Hour = AA{1}(6:end);
        Sec = AA{2}(5:end);
        Action = AA{5}(9:end);
        if isequal(Action,'*ARMED**'); continue; end
        if isempty(EventLog.CapControl.(Name).Regulations(1).time)
            cid = 1;
        else
            cid = length(EventLog.CapControl.(Name).Regulations) + 1;
        end
        EventLog.CapControl.(Name).Regulations(cid).time = (str2double(Hour)*3600)+str2double(Sec);
        EventLog.CapControl.(Name).Regulations(cid).action = Action;
    end
end

if ~isempty(EventLog.RegControl)
    Names_Reg = fieldnames(EventLog.RegControl);
    A = EventLog.RegControl;
    Reg = zeros(Pts,length(Names_Reg));
    isReg = 1;
else
    isReg = 0;
end
if ~isempty(EventLog.CapControl)
    Names_Cap = fieldnames(EventLog.CapControl);
    B = EventLog.CapControl;
    Cap = zeros(Pts,length(Names_Cap));
    isCapc = 1;
else
    isCapc = 0;
end

%% Treat Regcontrol operations
if isReg==1
    for i=1:length(Names_Reg)
        RegCtr = Names_Reg{i};
        %% tap position
        act = regexp({EventLog.RegControl.(RegCtr).Regulations.action}','\ ','split');
        x = zeros(length(act),2);
        t = [EventLog.RegControl.(RegCtr).Regulations.time]';
        for j = 1:length(act)
            x(j,:) = [str2double(act{j}{2}) str2double(act{j}{5})];
        end
        % voltage change per tap change
        if size(x,1)>1
            dTap = (x(2,2)-x(1,2))/x(2,1);
        else
            dTap = 0.0625;
        end
        tapPos = [x zeros(length(act),1)];
        % add initial tap position if not specified
        if EventLog.RegControl.(RegCtr).Regulations(1).time ~= 0
            % add intial time to time series
            t = [0;t];
            % add initial voltage value
            tapPos = [0 tapPos(1,2)-dTap*tapPos(1,1) 1; tapPos];
            % calculate tap position in integer (-15 to 16 for example)
            for k = 2:size(tapPos,1)
                tapPos(k,3) = tapPos(k-1,3)+tapPos(k,1);
            end
        else
            % TODO: this might need some more work
            
        end
        clear T V;
        %     Reg(1,i)=TotalA-1;
        tap_position = struct; T = []; V = [];
        for j=1:length(A.(RegCtr).Regulations)
            indexA = (A.(RegCtr).Regulations(j).time/step)+1;
            %counds all regcontrol events. if several events at one time step it adds 1 event to the
            %existing number at that time steps. it might take OpenDSS several regcontrol events for each control iteration
            %to find best solution.
            %         Reg(indexA,i) = Reg(indexA,i)+1;
            %counds one regcontrol event per time step. even if OpenDSS requires more
            %contol iterations to solve a step only one event is counted. in the real
            %world the voltage regulator would just jump several taps at a particular
            %time rather than go back and forth looking for best solution. this aproach
            %is more realistic than the one above
            Reg(indexA,i) = 1;
            tap_position(1,j).(RegCtr) = (EventLog.RegControl.(RegCtr).Regulations(j).time);
            T(j)=tap_position(1,j).(RegCtr)/30;
            tap_position(2,j).(RegCtr) = str2double(act{j}{5});
            V(j)=tap_position(2,j).(RegCtr);
        end
        % clean data
        timeCleaned = unique(T);
        clear ValueCleaned vartmp;
        for diffTi=1:length(timeCleaned)
            x = find(ismember(T,timeCleaned(diffTi)));
            ValueCleaned(diffTi)=V(x(end));
        end
        
        AB = zeros(1,Pts);
        for a=1:length(timeCleaned)-1
            AB(timeCleaned(a):Pts) = zeros(1,Pts-timeCleaned(a)+1)+ValueCleaned(a);
        end
        tap2plot.(RegCtr) = AB;
        try
            Cap(1,i) = sum(Cap(2:end,i));
            
            Reg(1,i) = sum(Reg(2:end,i));
        end
    end
    
    %sum all events by hour and add to the end of Reg matrix
    Reg(2:end,end)=sum(Reg(2:end,:)')';
else
    Reg='';
    tap2plot='';
end

%% Treat Capcontrol operations
if isCapc==1
    
    for i=1:length(Names_Cap)
        CapCtr = Names_Cap{i};
        act = regexp({B.(CapCtr).Regulations.action}','\ ','split');
        clear T V;
        T(1)=1;
        V(1)=0;
        Cap(1,i)=TotalB-1;
        cap_position(1,1).(CapCtr)= 1;
        cap_position(2,1).(CapCtr) = 0;
        for j=1:length(B.(CapCtr).Regulations)
            indexB = (RegulationsB(j).time/step)+1;
            %         Cap(indexB,i) = Cap(indexB,i)+1;
            %see comment above regarding countin all regcontrol events. same applies to
            %capcontol events.
            
            Cap(indexB,i) = 1;
            cap_position(1,j).(CapCtr) = (EventLog.CapControl.(CapCtr).Regulations(j).time);
            T(j)=cap_position(1,j).(CapCtr)/30;
            switch act{j}{1}
                case '*STEP'
                    switch act{j}{2}
                        case 'UP**'
                            V(j)=V(j-1) +1;
                        case 'DOWN**'
                            V(j)=V(j-1) -1;
                    end
                    
                case'*OPENED**'
                    V(j)= V(j-1)-1;
                    NumCapStep = V(j);
                    
                case '*CLOSED**'
                    V(j)=V(j-1) +1;
            end
            
        end
        V=V+circuit.capacitor(ismember(lower({circuit.capacitor.Name}),Names_Cap{i})).Numsteps;
        
        % clean data to create an array of value
        timeCleaned = unique(T);
        clear ValueCleaned vartmp;
        for diffTi=1:length(timeCleaned)
            act = find(ismember(T,timeCleaned(diffTi)));
            ValueCleaned(diffTi)=V(act(end));
        end
        
        AB = zeros(1,Pts);
        for a=1:length(timeCleaned)-1
            AB(timeCleaned(a):Pts) = zeros(1,Pts-timeCleaned(a)+1)+ValueCleaned(a);
        end
        cap2plot.(CapCtr) = AB;
        Cap(1,i) = sum(Cap(2:end,i));
    end
    %sum all events by hour and add to the end of Cap matrix
    Cap(2:end,end)=sum(Cap(2:end,:)')';
else
    Cap = '';
    cap2plot = '';
end
if ~exist('Reg', 'var'); Reg=zeros(1,2880);end
if ~exist('Cap', 'var'); Cap=zeros(1,2880);end
end

