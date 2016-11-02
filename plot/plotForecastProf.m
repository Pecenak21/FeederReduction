function plotForecastProf(forecast, dayid, fName, fcProfileId, outDir)
% plot all forecast profiles
if exist('forecast','var') 
    fc = forecast;
else
    fc = loadForecast(dayid, fName, fcProfileId);
end
if exist('fcProfileId','var') && ~isempty(fcProfileId)
    id = fcProfileId;
else
    id = 1:size(fc.profile,2);
end

if ~exist('outDir','var') || isempty(outDir)
    outDir = 'tmp';
end
if ~exist(outDir,'dir'), mkdir(outDir); end

for i = 1:length(id)
    h = figure(30); plot(fc.time,fc.profile(:,id(i))); datetick;
    title(['Site # ' num2str(id(i))]);
    fn = [outDir '/fig_forecast_' fName '_' dayid '_'  num2str(id(i)) '.png'];
    saveFigure(h,fn);
end

end