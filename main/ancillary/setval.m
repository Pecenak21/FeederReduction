function setval(dssEngine,component,field,valInString)
% examples: 
%           setval(o,c.transformer(1),'taps','[1 1.1]')
%           setval(o,'Transformer.520_838G','taps','[1 1.1]')
if ~isempty(strfind(class(component),'dss'))
    cName = [class(component) '.' component.Name];
    cName = cName(4:end);
elseif ~isempty(strfind(class(component),'Interface.OpenDSS_Engine'))
    type = class(component); 
    type = strsplit(type,'.'); type = type{end};
    type = type(2:end-1);
    cName = [type '.' component.Name];
elseif ischar(component)
    cName = component;
elseif isempty(component) % assume setting some special value for the engine
    cName = '';
else
    error('Component input is invalid. See doc (''doc getval'') for usage examples.');
end

if ~isempty(cName), 
    cmd = [cName '.' field '=' valInString];
else
    cmd = ['set ' field '=' valInString];
end
dssEngine.Text.Command = cmd;
if ~isempty(dssEngine.Text.Result)
    error(dssEngine.Text.Result)
end

end