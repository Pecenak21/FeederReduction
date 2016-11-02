function c = addExistingPV(c,feederId,data)
% add existing pvs to the feeder. Only supports 5 SDG&E feeder for now.

if ischar(feederId), feederId = lower(feederId); end
switch feederId
    case {'pointloma','480','cabrillo','b',480}
        pv = data.Point_Loma;
    case {'fallbrook','520','avocado','a',520}
        pv = data.Fallbrook;
    case {'alpine','355','d',355}
        pv = data.Alpine;
    case {'valleycenter','909','c',909}
        pv = data.Valley_Center;
    case {'ramona','971','creelman','e',971}
        pv = data.Ramona;
    otherwise
        error('Not supported feeder. Please check or write new code to support this feeder.');
end

% clean up data, removing invalid entries
pv = pv(~isnan([pv.lat]));

% Classify by Group B and A&C (2 sets of pv data from different sources)
% group AC's address is real address obtained from CSI data
% group B's address is in lat,lon getting from SDG&E 1-line PV maps
groupAC = [strfind([pv.Group],'A') strfind([pv.Group],'C')];
groupB = strfind([pv.Group],'B');

% combine sites if they are at the same location
[A, B]=ismember({pv(groupAC).Address}',{pv(groupB).Address}');
% Change sq m & kw datas
for i=1:length(groupAC)
    if A(i)==1
        pv(i).Size_kW_ = pv(B(i)+length(groupAC)).Size_kW_;
        pv(i).Area_sq_m_ = pv(B(i)+length(groupAC)).Area_sq_m_;
    end
    if A(i)==0
        pv(i).Size_kW_ = pv(i).DC;
    end
end
% Remove the replicated sites after combining
pv(length(groupAC)+B(B>0)) = [];

%% Auto find bus for each of the PV site based on its location
p=dsspvsystem;
for i=1:length(pv)
    % locate bus for the PV system by finding closest bus to it
    busId = lower(locateBus(c,pv(i).lon,pv(i).lat)); 
    
    % construct the pvsystem
    p(i).Name = [num2str(i) 'PV_' busId];
    p(i).kVA = pv(i).Size_kW_;
    p(i).Pmpp = pv(i).Size_kW_;
    [p(i).bus1, p(i).phases] = fixBus(c,busId);
end

c.pvsystem = p;
c = convertPVModel(c); % convert to model using power as input and output

end