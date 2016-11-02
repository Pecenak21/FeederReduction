function c = mvCompCloser(c,comp,numStep)
% move components closer to main line based on bus connections.
% only tested numStep=1 currently. For numStep > 1, please do check results.
%
% input
%       c:        circuit
%       comp:     component name
%       numStep:  number of steps to get closer to the main line
%
% example:
%       c = mvCompCloser(c,'pvsystem',1);

if ~exist('numStep','var') || isempty(numStep)
    numStep = 1;
end
for j = 1:numStep
    for i = 1:length(c.(comp))
        [~,lid] = ismember(c.(comp)(i).bus1,{c.line.bus2});
        if isempty(lid)
            disp(i);
        end
        c.(comp)(i).bus1 = c.line(lid).bus1;
    end
    % move one more step closer to the main line
    for i = 1:length(c.(comp))
        [~,lid] = ismember(c.(comp)(i).bus1,{c.line.bus2});
        if ~isempty(lid)
            c.(comp)(i).bus1 = c.line(lid).bus1;
        else
            disp(i);
        end
    end
end

end