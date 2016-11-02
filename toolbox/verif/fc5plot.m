function fc5plot(USI_ID, date, site_ID)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('Warning: Make sure your folder contains consecutive files!')

% Load UCSD Deployment
javaCheckInit();
target = siDeployment('UCSD');

% If no site_ID specified, use average DEMROES data
if nargin == 2
	% Average data
	site_ID = 1;
	site_end = 6;
	site_name = 'Average';
	fprintf('===%s Data Mode===\n', site_name)
elseif nargin == 3
	% Site specific data
	site_end = site_ID;
	site_name = target.footprint.GHInames{site_ID};
	fprintf('===Site-specific Data Mode for %s===\n', site_name)
else
	error('Incorrect input: ncplot(USI ID, yyyymmdd, DEMROES site ID)')
end

% Force date to be a string for sprintf purposes
if ~ischar(date)
	date = num2str(date);
end

% Assign directories
inputDir = sprintf('C:/Users/Handa/Desktop/ForecastOut/USI_1_%i/%s/',USI_ID,date);
outputDir = 'C:/Users/Handa/Documents/Research/NowcastPlots/';
if(exist(outputDir,'dir') == 0)
	mkdir(outputDir);
end

% Load files
% Forecast files
forecast_files = dir(sprintf('%s/forecast/', inputDir));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
	% Average data
	fprintf('Loading %s data forecasts...\n', site_name)
	
	% Initialize
	first_forecast = load(sprintf('%s/forecast/%s', inputDir, forecast_files(3).name));
	timeUTC(1) = first_forecast.time(31);
	forecast_kt(1) = first_forecast.ktavg(31);
	fprintf('Start time: %s UTC\n', datestr(timeUTC(1)))
	
	% Create UTC time vector
	for z = 2:numel(forecast_files)-2
		timeUTC(z) = timeUTC(z-1) + 30*1/24/60/60;
	end
	
	for t = 4:numel(forecast_files)
% 		fprintf('Loading %ith file of %i...\n', t, numel(forecast_files)-1)
		forecast = load(sprintf('%s/forecast/%s', inputDir, forecast_files(t).name), 'ktavg', 'time');
		forecast_kt(t-2) = forecast.ktavg(31);
	end
elseif nargin == 3
	% Site-specific data
	fprintf('Loading forecasts for %s...\n', site_name)
	
	% Initialize
	first_forecast = load(sprintf('%s/forecast/%s', inputDir, forecast_files(3).name));
	timeUTC(1) = first_forecast.time(31);
	forecast_kt(1) = first_forecast.station(site_ID).ktavg(31);
	fprintf('Start time: %s UTC\n', datestr(timeUTC(1)))
	
	% Create UTC time vector
	for z = 2:numel(forecast_files)-2
		timeUTC(z) = timeUTC(z-1) + 30*1/24/60/60;
	end
	
	for t = 4:numel(forecast_files)
% 		fprintf('Loading %ith file of %i...\n', t, numel(forecast_files))
		forecast = load(sprintf('%s/forecast/%s', inputDir, forecast_files(t).name), 'station');
		forecast_kt(t-2) = forecast.station(site_ID).ktavg(31);
	end
end

fprintf('End time: %s UTC\n', datestr(timeUTC(numel(timeUTC))))

% Create PST time vector

for n = 1:z
	timePST(n) = addtodate(timeUTC(n), -8, 'hour');
end

% Get day of year
StartDate = datenum(first_forecast.time(1));
PrevYear = datenum(str2double(datestr(StartDate, 'yyyy'))-1, 12,31);
DOY = floor(StartDate-PrevYear);

% Preallocate GHI matrix
GHI(1:numel(timePST),1:6) = NaN;

% Load GHI data
for ID = site_ID:site_end
	tic
	fprintf('Loading GHI data...\n')
	GHIdata(ID) = load(siNormalizePath(sprintf('%s/%s/%s_%i_%i.mat','$KleisslLab4TB1/database/DEMROES/1s_by_DoY', target.footprint.GHInames{ID}, target.footprint.GHInames{ID}, str2double(datestr(first_forecast.time(1), 'yyyy')), DOY)));
	GHIstart = find(abs(GHIdata(ID).time_day(:,1) - timePST(1)) < datenum([0000 00 00 00 00 01])); % Search for start time in DEMROES data and convert into UTC time. DEMROES -8 hours w.r.t. UTC.
	GHIend = find(abs(GHIdata(ID).time_day(:,1) - timePST(numel(timePST))) < datenum([0000 00 00 00 00 01]));
	% Concatenate
	GHI(:,ID) = GHIdata(ID).GHI_day(GHIstart:30:GHIend);
	GHITIME(:,ID) = GHIdata(ID).time_day(GHIstart:30:GHIend);
	toc
end

GHIavg = zeros(length(timeUTC),1);
for n = 1:length(GHI(:,1))
	GHIavg(n) = nanmean(GHI(n,:));
end

% Convert DEMROES GHI Data to kt
csk = clearSkyIrradiance( target.ground.position , timeUTC, target.tilt, target.azimuth );
GHI_kt = GHIavg./csk.gi;

% Plot nowcast kt onto DEMROES kt
a = figure;
set(a, 'Position', [50 50 1500 900])
hold on
plot(timeUTC, forecast_kt, 'k')
plot(timeUTC, GHI_kt, 'g')
datetick('x')
axis tight
xlabel('Time (UTC) [hh:mm]')
ylabel('Clear Sky Index kt [-]')
legend('Forecast kt', sprintf('%s DEMROES calculated kt', site_name), 'Location', 'Best')
title(sprintf('5 minute forecast kt overlaid on DEMROES calculated kt for %s', datestr(timeUTC(1), 'mmm-dd-yyyy')))

end