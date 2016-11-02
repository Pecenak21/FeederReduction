function c = addEnergyStorage(c,comp)
% add ES to all nodes that the specified components are connected to
%
% Note: new ES object will have
%       kWhrated = comp's kW
%       kWrated = comp's kW*.2
%
% input
%       c:      circuit
%       comp:   components in circuit (i.e. 'pvsystem', 'generator')
%
% example
%       c = addEnergyStorage(c,c.pvsystem);

if strfind(class(comp),'dss')
    g = comp; type = strrep(class(comp),'dss','');
elseif ischar(comp)
    g = c.(comp); type = lower(g);
else
    error('not supported ''comp'' input!');
end

es(length(g)) = dssstorage;
for i = 1:length(g)
    es(i).Name = sprintf('es_%d',i);
    es(i).bus1 = g(i).bus1;
    es(i).phases = g(i).phases;
    switch(type)
        case 'generator'
            es(i).kWhrated = g(i).kW;
            es(i).kWrated = g(i).kW*.2;
        case 'pvsystem'
            es(i).kWhrated = g(i).kVA/1.1;
            es(i).kWrated = es(i).kWrated*.2;
        otherwise
            error('not supported ''comp'' yet. please write code to support it here');
    end
    es(i).reserve = 0;
    es(i).EffCharge = 100;
    es(i).EffDischarge = 100;
    es(i).Idlingkvar = 0;
    es(i).IdlingkW = 0;
    
    % initial set up: 100% full and dischargin
    es(i).stored = 100; % default 100%
    es(i).State = 'idling';
end

c.storage = es;
end