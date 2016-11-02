function [ val ] = dssgetelem( dssengine, element )
%DSSGETELEM return element of circuit after simulation
%
% dssengine: dssEngine after simulation or first output of dssget  
% element: element(s) of interest either in dssObj or in string format: TYPE.NAME 
%         exp: c.fault(1) or 'fault.1'

o = dssengine;

if ~ischar(element) && strfind(class(element),'dss')
    c = class(element);
    element = [c(4:end) '.' element.Name];
end

o.text.command = ['select ' element];
val = dssengine.ActiveCircuit.ActiveCktElement;

end

