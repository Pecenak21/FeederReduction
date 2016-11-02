function [EventLog Reg Cap] = EventLogEvaluation_Vadim(Event, Pts, step)

%evaluated the EventLog for the current OpenDSS simulation. 
%INPUTS: 
%"Event" - event log (capacitor switching, regulator tap changes)  
%"Pts" - number of 'solves' (e.g. 24 for a daily
%simulaiton with 1h stepsize)

%OUTPUTS:
%"EventLog" - structure containing all the events sorted by control devices,
%hour of simulaiton, and performed action
%"Reg" - array summarizing RegControl events. first row represents the total
%events for the corresponding controller. The remaining rows represent the 
%point of time in the current simulation. The columns represent the 
%individual controllers lined up in the following manner:
% 'I' 'II' 'III' 'IV' 'V' 'VI' 'Sub'.
%"Cap" - array summarizing CapControl events (same structure as "Reg").

disp('- Creating EventLog');
EventLog = struct('RegControl', {''}, 'CapControl', {''});

EventLog.RegControl.reg1 = struct('ID', 'Reg_I', 'Regulations', {''});
EventLog.RegControl.reg1.Regulations = struct('time', 0, 'action', 'Nada');

% EventLog.RegControl.t_804 = struct('ID', 'Reg_II', 'Regulations', {''});
% EventLog.RegControl.t_804.Regulations = struct('time', 0, 'action', 'Nada');
% 
% EventLog.RegControl.t_805 = struct('ID', 'Reg_III', 'Regulations', {''});
% EventLog.RegControl.t_805.Regulations = struct('time', 0, 'action', 'Nada');
% 
% EventLog.RegControl.t_811 = struct('ID', 'Reg_IV', 'Regulations', {''});
% EventLog.RegControl.t_811.Regulations = struct('time', 0, 'action', 'Nada');
% 
% EventLog.RegControl.t_707 = struct('ID', 'Reg_V', 'Regulations', {''});
% EventLog.RegControl.t_707.Regulations = struct('time', 0, 'action', 'Nada');
% 
% EventLog.RegControl.t_708 = struct('ID', 'Reg_VI', 'Regulations', {''});
% EventLog.RegControl.t_708.Regulations = struct('time', 0, 'action', 'Nada');
% 
% EventLog.RegControl.t_712 = struct('ID', 'Reg_VII', 'Regulations', {''});
% EventLog.RegControl.t_712.Regulations = struct('time', 0, 'action', 'Nada');
% 
% EventLog.RegControl.t_729 = struct('ID', 'Reg_VIII', 'Regulations', {''});
% EventLog.RegControl.t_729.Regulations = struct('time', 0, 'action', 'Nada');
% 
% EventLog.RegControl.t_730 = struct('ID', 'Reg_IX', 'Regulations', {''});
% EventLog.RegControl.t_730.Regulations = struct('time', 0, 'action', 'Nada');
% 
% EventLog.RegControl.t_755 = struct('ID', 'Reg_X', 'Regulations', {''});
% EventLog.RegControl.t_755.Regulations = struct('time', 0, 'action', 'Nada');
% 
% EventLog.RegControl.t_758 = struct('ID', 'Reg_XI', 'Regulations', {''});
% EventLog.RegControl.t_758.Regulations = struct('time', 0, 'action', 'Nada');
% 
% EventLog.RegControl.t_759 = struct('ID', 'Reg_XII', 'Regulations', {''});
% EventLog.RegControl.t_759.Regulations = struct('time', 0, 'action', 'Nada');
% 
% EventLog.RegControl.t_765 = struct('ID', 'Reg_XIII', 'Regulations', {''});
% EventLog.RegControl.t_765.Regulations = struct('time', 0, 'action', 'Nada');

EventLog.RegControl.san_berdo = struct('ID', 'Reg_Sub', 'Regulations', {''});
EventLog.RegControl.san_berdo.Regulations = struct('time', 0, 'action', 'Nada');

EventLog.CapControl.cap_4250412e = struct('ID', 'Cap_I', 'Regulations', {''});
EventLog.CapControl.cap_4250412e.Regulations = struct('time', 0, 'action', 'Nada');

for i=1:length(Event)
    AA = strread(Event{i},'%s','delimiter',',');
    if strcmpi(AA{4}(9),'R')
        Name = AA{4}(19:end);
        if ~isfield(EventLog.RegControl, Name)
            disp(['The RegControl ' Name ' is not on the current list. Add it NOW']);
            return
        end
        Hour = AA{1}(6:end);
        Sec = AA{2}(5:end);
        Action = AA{5}(9:end-1);
        EventLog.RegControl.(Name).Regulations(end+1).time = (str2num(Hour)*3600)+str2num(Sec);
        EventLog.RegControl.(Name).Regulations(end).action = Action;
        
    else
        Name = AA{4}(19:end);
        if ~isfield(EventLog.CapControl, Name)
            disp(['The CapControl ' Name ' is not on the current list. Add it NOW']);
            return
        end
        Hour = AA{1}(6:end);
        Sec = AA{2}(5:end);
        Action = AA{5}(9:end);
        EventLog.CapControl.(Name).Regulations(end+1).time = (str2num(Hour)*3600)+str2num(Sec);
        EventLog.CapControl.(Name).Regulations(end).action = Action;
    end
end

Names_Reg = fieldnames(EventLog.RegControl);
Names_Cap = fieldnames(EventLog.CapControl);
A = EventLog.RegControl;
B = EventLog.CapControl;
Reg = zeros(Pts+1,length(Names_Reg)+1);
Cap = zeros(Pts+1,length(Names_Cap)+1);
for i=1:length(Names_Reg)
    RegCtr = Names_Reg{i};
    RegulationsA = A.(RegCtr).Regulations;
    TotalA = length(RegulationsA);
%     Reg(1,i)=TotalA-1;
    for j=2:TotalA
        indexA = (RegulationsA(j).time/step)+1;
%counts all regcontrol events. if several events at one time step it adds 1 event to the
%existing number at that time steps. it might take OpenDSS several regcontrol events for each control iteration
%to find best solution.
%         Reg(indexA,i) = Reg(indexA,i)+1;
%counts one regcontrol event per time step. even if OpenDSS requires more
%contol iterations to solve a step only one event is counted. in the real
%world the voltage regulator would just jump several taps at a particular
%time rather than go back and forth looking for best solution. this aproach
%is more realistic than the one above
        Reg(indexA,i) = 1;
    end
    Reg(1,i) = sum(Reg(2:end,i));
end
%sum all events by hour and add to the end of Reg matrix
Reg(2:end,end)=sum(Reg(2:end,:)')';

for i=1:length(Names_Cap)
    CapCtr = Names_Cap{i};
    RegulationsB = B.(CapCtr).Regulations;
    TotalB = length(RegulationsB);
%     Cap(1,i)=TotalB-1;
    for j=2:TotalB
        indexB = (RegulationsB(j).time/step)+1;
%         Cap(indexB,i) = Cap(indexB,i)+1;
%see comment above regarding counting all regcontrol events. same applies to
%capcontol events.
        Cap(indexB,i) = 1;
    end
    Cap(1,i) = sum(Cap(2:end,i));
end
%sum all events by hour and add to the end of Cap matrix
Cap(2:end,end)=sum(Cap(2:end,:)')';

end