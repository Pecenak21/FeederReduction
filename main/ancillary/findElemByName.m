function [elem, id] = findElemByName(c,compType,name)
name = lower(name);
id = find(ismember(lower([c.(compType).Name]),name));
elem = [];
if ~isempty(id)
    elem = c.(compType)(id);
end
end