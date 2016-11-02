%% check forecast profile
dayid = '20121219'; fName = 'fallbrook'; fcProfileId = []; outDir = 'Y:/database/gridIntegration/PVimpactPaperFinal/figForecast';
plotForecastProf(dayid, fName, fcProfileId, outDir);
%%
dayid = '20121218'; fName = 'fallbrook'; fcProfileId = []; 
outDir = ['Y:/database/gridIntegration/PVimpactPaperFinal/'...
            'figForecast/' fName '_' dayid];
plotForecastProf(dayid, fName, fcProfileId, outDir);
%% plot all forecast profiles
dayid = {'20121219','20121218','20121214'};
fName = {'fallbrook','pointloma','valleycenter','ramona','alpine'}; fcProfileId = [];
simDir = 'Y:/database/gridIntegration/PVimpactPaperFinal/';

for i = 1:length(dayid)
    day = dayid{i};
    for j = 1:length(fName)
        n = fName{j};
        outDir = [simDir 'figForecast/' n '_' day];
        plotForecastProf([], day, n, fcProfileId, outDir);
    end
end

%% generate new forecast profile for all feeders for all days
dayid = {'20121219','20121218','20121214'};
deploySite = {'Fallbrook_432pvs_UCSD'...
    'PointLoma_340pvs_UCSD'...
    'ValleyCenter_387pvs_UCSD'...
    'Ramona_387pvs_UCSD'...
    'Alpine_364pvs_UCSD'};
fName = {'fallbrook','pointloma','valleycenter','ramona','alpine'}; 

for i = 1:length(dayid)
    for j = 1:length(deploySite)
        fc = loadForecast(dayid{i}, fName{j});
        fillForecastProfile(fc, deploySite{j}, fName{j},dayid{i});
    end
end

%% plot all filled forecast profiles
dayid = {'20121219','20121218','20121214'};
fName = {'fallbrook','pointloma','valleycenter','ramona','alpine'}; 
deploySite = {'Fallbrook_432pvs_UCSD'...
    'PointLoma_340pvs_UCSD'...
    'ValleyCenter_387pvs_UCSD'...
    'Ramona_387pvs_UCSD'...
    'Alpine_364pvs_UCSD'};
simDir = 'Y:/database/gridIntegration/PVimpactPaperFinal/';
fcProfileId = [];
filled = 1; % use filled forecast profiles
for i = 1:length(dayid)
    for j = 1:length(fName)
        outDir = [simDir 'figForecastFilled/' n '_' day];
        fc = fillForecastProfile(fc, deploySite{j}, fName{j},dayid{i});
        plotForecastProf(fc, dayid{i}, fName{j}, fcProfileId, outDir);
    end
end