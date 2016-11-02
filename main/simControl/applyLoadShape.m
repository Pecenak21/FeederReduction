function c = applyLoadShape(c,timeId,comp,loadshapeId)
global indent; if isempty(indent), indent = ''; end

if strfind(class(comp),'dss')
    g = comp; type = strrep(class(comp),'dss','');
elseif ischar(comp)
    g = c.(comp); type = lower(g);
else
    error('not supported ''comp'' input!');
end

switch type
    case 'pvsystem'
        fn = 'pmpp';
    case 'storage'
end
for i = 1:length(c.comp)
    c.comp(i).fn = c.loadshape(loadshapeId(i)).Mult(timeId);
end
end