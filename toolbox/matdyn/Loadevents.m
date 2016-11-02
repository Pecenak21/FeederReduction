function [event,buschange,linechange]=Loadevents(casefile_ev)
      
% [event,buschange,linechange]=Loadevents(casefile_ev)
% 
% Loads event data
% 
% INPUTS
% casefile_ev = m-file with event data
% 
% OUTPUTS
% event = list with time and type of events
% buschange = time, bus number, bus parameter, new value
% linechange = time, branch number, branch parameter, new value


% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%%
if isstruct(casefile_ev)
    event = casefile_ev.event;
    type1 = casefile_ev.buschange;
    type2 = casefile_ev.linechange;
else
    [event,type1,type2] = feval(casefile_ev);
end


buschange=zeros(size(event,1),4);
linechange=zeros(size(event,1),4);

i1=1;
i2=1;

for i=1:size(event,1)
    if event(i,2)==1
        buschange(i,:) = type1(i1,:);
        i1=i1+1;     
	elseif event(i,2)==2
        linechange(i,:) = type2(i2,:);
        i2=i2+1;
    end
end
        
return;