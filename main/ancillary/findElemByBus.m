function [elem, elemId] = findElemByBus(c,bus,careAboutPhase,compType)
% find element(s)/ component(s) in the circuit that connect(s) to specified bus
%
% input:
%       c:          circuit struct or struct array of components/ elems
%       bus
%       compType:   only needed when circuit is input
%       careAboutPhase:  if care about phase then only bus1.1.2 matches
%                       bus1.1.2, else bus1.1 matches bus1.1.2.3 (phase-careless)

if isa(c,'COM.OpendssEngine_dss')
    n = c.activeCircuit.AllElementNames;
    elem = {}; elemId = [];
    for i = 1:length(n)
        if ismember(lower(bus),lower(cleanBus(getval(c,n{i},'busnames'))))
            elem = [elem,n{i}];
            elemId = [elemId, i]
        end
    end
return;
end
if ~exist('careAboutPhase','var') || isempty(careAboutPhase)
    careAboutPhase = 0;
end
elem = []; cid = [];
if isstruct(c) && isfield(c,'circuit')
    compNames = fieldnames(c);
    
    for comId = 1:length(compNames)
        cName = compNames{comId};
        if strfind(class(c.(cName)),'dss')
            [x, xid] = findElem(c.(cName),bus);
            if ~isempty(x)
                elem.(cName) = x;
                elem.([cName '_id']) = xid;
            end
        end
    end
elseif strfind(class(c),'dss')
    [elem, elemId] = findElem(c,bus);
elseif exist('compType','var') && isfield(c,compType)
    [elem, elemId] = findElem(c.(compType),bus);
end

    function [clist, cId] = findElem(component,bus)
        bus = lower(bus);
        cId = [];
        bStrList = {'bus','bus1','bus2','buses'};
        fn = lower(fieldnames(component));
        bStr = ismember(bStrList,fn);
        bStr = bStrList(bStr);
        if ~isempty(bStr)
            for i = 1:length(bStr)
                if ismember(bStr{i},{'bus','buses'}) && isempty(strfind(class(component),'regcontrol'))
                    if strcmp(bStr{i},'bus')
                        if isempty([component.bus]), continue; end
                        blist = lower(reshape([component.bus],2,length(component))');
                    else
                        if isempty([component.buses]), continue; end
                        blist = lower(reshape([component.buses],2,length(component))');
                    end
                    if ~careAboutPhase
                        blist(:,1) = cleanBus(blist(:,1));
                        blist(:,2) = cleanBus(blist(:,2));
                    end
                    y = [find(ismember(blist(:,1),bus)) find(ismember(blist(:,2),bus))];
                else
                    blist = lower([component.(bStr{i})]);
                    if ~careAboutPhase
                        blist = cleanBus(blist);
                    end
                    if ischar(blist), blist = {blist}; end
                    y = find(ismember(blist,bus));
                end
                cId = [cId y];
            end
            cId = unique(cId);
            clist = component(cId);
        else
            clist = [];
        end
    end

end