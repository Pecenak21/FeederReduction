function [c, mult] = applyPenLevel(c, pen, excludeScalingUpIds, currentPen)
% change penetration level
% currently only supports PVsystem penetration
%
% input: 
%           excludeScalingUpIds     : systems' ids not to scale up (e.g. 2 big PV systems at Fallbrook)

if exist('excludeScalingUpIds','var') && ~isempty(excludeScalingUpIds)
    excludeId = excludeScalingUpIds;
else
    excludeId = [];
end
includeId = 1:length(c.pvsystem);
type = 'pvsystem';
% calculate current penetration level (total pv peak output in KW/ total load in KW)
% assume all systems are enabled

if strcmpi(type,'pvsystem')
    if exist('currentPen','var')
        cPen = currentPen;
    else
        cPen = sum([c.pvsystem.Pmpp])/ sum([c.load.kw])*100;
    end
    
    % multiplication factor
    mult = pen/cPen;
    
    % if excluding some systems and only when scaling them up (mult > 1)
    if mult > 1 && ~isempty(excludeId)
        includeId = setdiff(1:length(c.pvsystem),excludeId);
        mult = (pen*sum([c.load.kw])/100 - sum(cell2mat([c.pvsystem(excludeId).Pmpp])) )...
            /sum(cell2mat([c.pvsystem(includeId).Pmpp]));
    end
    
    for i = 1:length(includeId)
        id = includeId(i);
        switch lower(type)
            case 'pvsystem'
                c.(type)(id).Pmpp = mult * c.(type)(id).Pmpp;
                c.(type)(id).kVA = c.(type)(id).Pmpp*1.05;
            otherwise
                error('Doesn''t handle other types besides pvsystem yet!');
        end
    end
end

end