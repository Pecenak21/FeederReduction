function [ o ] = setSolutionTime( o, time, startTime )
global oStartTime; 
% calculate solution time in seconds
if ~exist('startTime','var') || isemtpy(startTime)
    if ~isempty(oStartTime)
        startTime = oStartTime;
    else
        startTime = 0;
    end
end
time = time - startTime;
o.ActiveCircuit.Solution.dblHour = time*24;
end

