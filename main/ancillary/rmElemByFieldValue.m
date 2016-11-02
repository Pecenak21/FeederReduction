function c = rmElemByFieldValue(c,compType,field,val)
id = find(ismember([c.(compType).(field)],val));
if ~isempty(id)
    disp('Removing below components:');
    char(c.(compType)(id))
    c.(compType)(id) = [];
end
end