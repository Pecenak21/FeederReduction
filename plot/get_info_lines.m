function [pv,bus1,bus2,substation,LinesToPlot,f] = get_info_lines(feeder)
%% Declaring variables to be used
LinesToPlot = feeder.line;
busIds  = lower(feeder.buslist.id);
linesBus1 = lower(regexprep({LinesToPlot.bus1},'(\.[0-9]+)','')'); %take out the phase numbers on buses if they have them
linesBus2 = lower(regexprep({LinesToPlot.bus2},'(\.[0-9]+)','')'); %take out the phase numbers on buses if they have them

%% Getting the real buses coordinates.
[busesRealData f] = realCoordAccurate(feeder);
allBusesCoordinates = zeros(length(busesRealData),2);
for jj=1:length(busesRealData)
    allBusesCoordinates(jj,:) = busesRealData(jj).coord;
end

%% Getting buses coordinates, for future plot
[~,indexBus1] = ismember(linesBus1,busIds);
bus1.coord  = allBusesCoordinates(indexBus1,:);
bus1.names = lower(feeder.buslist.id(indexBus1,:));
[~,indexBus2] = ismember(linesBus2,busIds);
bus2.coord = allBusesCoordinates(indexBus2,:);


%% Getting PV locations and power
pv.bus1= lower({feeder.pvsystem.bus1});
pv.kVA = [feeder.pvsystem.kVA];
[~,indexPv] = ismember(pv.bus1,busIds);
pv.coord = allBusesCoordinates(indexPv,:);
%% Substation
[~,indexSubs] = ismember('sourcebus',busIds);
substation = allBusesCoordinates(indexSubs,:);
end