function val = getval(dssEngine,component,field)
% examples: 
%           getval(o,c.transformer(1),'taps')
%           getval(o,c.transformer,'taps')
%           getval(o,'Transformer.520_838G','taps')

if ~isempty(strfind(class(component),'dss')) && length(component) > 1
    val = cell(length(component),1);
    for i = 1:length(component)
        val{i} = getvalElmt(dssEngine,component(i),field);
    end
elseif ~isempty(strfind(class(component),'Interface.OpenDSS_Engine'))
    type = class(component); 
    type = strsplit(type,'.'); type = type{end};
    type = type(2:end-1);
    component = [type '.' component.Name];
    val = getvalElmt(dssEngine,component,field);
else
    val = getvalElmt(dssEngine,component,field);
end
end

function val = getvalElmt(dssEngine,component,field)
if strfind(class(component),'dss')
    cName = [class(component) '.' component.Name];
    cName = cName(4:end);
elseif ischar(component)
    cName = component;
else
    error('Component input is invalid. See doc (''doc getval'') for usage examples.');
end

dssEngine.Text.Command = ['? ' cName '.' field];
val = dssEngine.Text.Result;

% if no such property exists then try using dssgetval to get it if avalable
if strcmpi(val,'property unknown')
    val = dssgetval(dssEngine,component,field);
end

end