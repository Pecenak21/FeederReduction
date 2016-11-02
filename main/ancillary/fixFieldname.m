function fixedFn = fixFieldname(fn,obj)
fns = fieldnames(obj);
id = find(ismember(lower(fns),lower(fn)));
if ~isempty(id)
    fixedFn = fns{id(1)};
else
    fixedFn = '';
end
end
