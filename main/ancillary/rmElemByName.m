function c = rmElemByName(c,compType,name)
id = find(ismember(lower([c.(compType).Name]),lower(name)));
if ~isempty(id)
    disp('Removing below components:');
    char(c.(compType)(id))
    c.(compType)(id) = [];
end
end