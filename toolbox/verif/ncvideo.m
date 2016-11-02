function ncvideo(USI_ID, date, site_ID)

%% Setup

warning('Warning: Make sure your folder contains consecutive files!')

% Load UCSD Deployment
javaCheckInit();
target = siDeployment('UCSD');

% Load config
conf = readConf(siGetConfPath('tools.conf'));

% If no site_ID specified, use average DEMROES data
if nargin == 2
	% Average data
	site_ID = 1;
	site_end = 6;
	site_name = 'Average';
	site_X = target.DEMROES_X;
	site_Y = target.DEMROES_Y;
	fprintf('===%s Data Mode===\n', site_name)
elseif nargin == 3
	% Site specific data
	site_end = site_ID;
	site_name = target.footprint.GHInames{site_ID};
	site_X = target.DEMROES_X(site_ID);
	site_Y = target.DEMROES_Y(site_ID);
	fprintf('===Site-specific Data Mode for %s===\n', site_name)
else
	error('Incorrect input: ncplot(USI ID, yyyymmdd, DEMROES site ID)')
end

% Force date to be a string for consistency
if ~ischar(date)
	date = num2str(date);
end

% Assign directories
inputDir = sprintf('%s/USI_1_%i/%s/', conf.inputDir, USI_ID, date);
outputDir = sprintf('%s/NowcastVideos/%s/', conf.outputDir, site_name);
if(exist(outputDir,'dir') == 0)
	mkdir(outputDir);
end

%% Load files

% Forecast files list
forecast_files = dir(sprintf('%s/forecast/', inputDir));

if nargin == 2
	% Average data
	fprintf('Loading %s data forecasts...\n', site_name)
	
	% Initialize
	first_forecast = load(sprintf('%s/forecast/%s', inputDir, forecast_files(3).name));
	timeUTC(1) = first_forecast.time(1);
	nowcast_kt(1) = first_forecast.ktavg(1);
	fprintf('Start time: %s UTC\n', datestr(timeUTC(1)))
	
	% Create UTC time vector
	for z = 2:numel(forecast_files)-2
% 		timeUTC(z) = timeUTC(z-1) + 30*1/24/60/60;
		% ^ Inaccurate, loses seconds over time
		% v Get time using forecast file name
		timeUTC(z) = datenum(forecast_files(z+2).name(10:23), 'yyyymmddHHMMSS');
	end
	
	% Commence mass loading
	for t = 4:numel(forecast_files)
% 		fprintf('Loading %ith file of %i...\n', t, numel(forecast_files)-1)
		forecast = load(sprintf('%s/forecast/%s', inputDir, forecast_files(t).name), 'ktavg', 'time');
		nowcast_kt(t-2) = forecast.ktavg(1);
	end
elseif nargin == 3
	% Site-specific data
	fprintf('Loading forecasts for %s...\n', site_name)
	
	% Initialize
	first_forecast = load(sprintf('%s/forecast/%s', inputDir, forecast_files(3).name));
	timeUTC(1) = first_forecast.time(1);
	nowcast_kt(1) = first_forecast.station(site_ID).ktavg(1);
	fprintf('Start time: %s UTC\n', datestr(timeUTC(1)))
	
	% Create UTC time vector
	for z = 2:numel(forecast_files)-2
% 		timeUTC(z) = timeUTC(z-1) + 30*1/24/60/60;
		% ^ Inaccurate, loses seconds over time
		% v Get time using forecast file name
		timeUTC(z) = datenum(forecast_files(z+2).name(10:23), 'yyyymmddHHMMSS');
	end
	
	% Commence mass loading
	for t = 4:numel(forecast_files)
% 		fprintf('Loading %ith file of %i...\n', t, numel(forecast_files))
		forecast = load(sprintf('%s/forecast/%s', inputDir, forecast_files(t).name), 'station');
		nowcast_kt(t-2) = forecast.station(site_ID).ktavg(1);
	end
end

fprintf('End time: %s UTC\n', datestr(timeUTC(numel(timeUTC))))

% Create PST time vector

for n = 1:z
	timePST(n) = addtodate(timeUTC(n), -8, 'hour');
end

%% Load GHI Data

% Get day of year
StartDate = datenum(first_forecast.time(1));
PrevYear = datevec(StartDate);
PrevYear = datenum(PrevYear(1)-1, 12,31);
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

% Average GHI data if needed
GHIavg = zeros(length(timeUTC),1);
for n = 1:length(GHI(:,1))
	GHIavg(n) = nanmean(GHI(n,:));
end

% Convert DEMROES GHI Data to kt
csk = clearSkyIrradiance( target.ground.position , timeUTC, target.tilt, target.azimuth );
GHI_kt = GHIavg./csk.gi;

%% Generate Video

% Create video object
vidObj = VideoWriter(sprintf('%s/%s_%s',outputDir, datestr(timeUTC(1), 'mmm-dd-yyyy'), site_name), 'MPEG-4');
vidObj.Quality = 90;
vidObj.FrameRate = 10;
open(vidObj);

% Initialize figure window
a = figure('Color', [1 1 1]);
set(a, 'Position', [50 50 1500 900])

% Plot large graph of nowcast kt overlaid on DEMROES kt
plot(timeUTC, nowcast_kt, 'k', 'LineWidth', 2)
hold on
plot(timeUTC, GHI_kt, 'g', 'LineWidth', 2)
hold off
set(gca, 'FontSize', 12)
datetick('x')
axis tight
xlabel('Time (UTC) [hh:mm]')
ylabel('Clear Sky Index kt [-]')
legend('Nowcast kt', sprintf('%s DEMROES kt', site_name), 'Location', 'SouthWest')
title(sprintf('Nowcast kt overlaid on DEMROES calculated kt for %s', datestr(timeUTC(1), 'mmm-dd-yyyy')))

% Initialize figure window
b = figure('Color', [1 1 1]);
set(b, 'Position', [50 50 1500 900])

% Create & assign color map
naturalskymap = [ 1 1 1; 0 .4 .8; .9 .9 .9; .7 .7 .7 ];

% Load last cloud motion file since it contains every cloud velocity
pcloudmotion = load(sprintf('%s/cloudmotion/cloudmotion_%s.mat', inputDir, datestr(timeUTC(end), 'yyyymmddHHMMSS')));

% Loop to incorporate video and zoomed graph
for pidx = 1:numel(timeUTC)
	% Load RBR and cloud decision from projection file
	% IMPORTANT: If projection file structure changes in the future, this
	% will need to be altered to accommodate
	ptime = datestr(timeUTC(pidx), 'yyyymmddHHMMSS');
	
	% Load pertinent forecast files
	pprojection = load(sprintf('%s/projection/projection_%s.mat', inputDir, ptime), 'rbr', 'dec');
	pforecast = load(sprintf('%s/forecast/forecast_%s.mat', inputDir, ptime), 'shadow');
	pcloudheight = load(sprintf('%s/cloudheight/cloudheight_%s.mat', inputDir, ptime));
		
	% Raw image (top left)
	subplot(2,3,1, 'replace');
	set(gca, 'Position', [0.015 0.5 0.3 0.47], 'FontSize', 12)
	img = imread(sprintf('%s/usi1-%i/%s/%s_prev.jpg', siNormalizePath('$KLEISSLLAB18-1/database/USI/images/'), USI_ID, date, ptime));
	imagesc(img)
	axis image
	axis off
	title('Raw Image')
	
	% RBR image (top middle)
	subplot(2,3,2, 'replace');
	set(gca, 'Position', [0.34 0.5 0.33 0.47], 'FontSize', 12)
 	mask = isnan(pprojection.rbr);
	h = imagesc(pprojection.rbr, [0.5 1.5]);
	colorbar
	axis image
	axis off
	title('Red-Blue-Ratio')
	
	% Cloud decision (top right)
	subplot(2,3,3, 'replace');
	set(gca, 'Position', [0.675 0.5 0.33 0.47], 'FontSize', 12)
	subimage(uint8(pprojection.dec+1), naturalskymap)
	axis image
	axis off
	title('Cloud Decision')
	
	% mask the RBR image
	c = get(h,'CData');
	c(repmat(mask,[1 1 3])) = 255;
	set(h,'CData',c);
		
	% Shadow map (bottom left)
	subplot(2,3,4)
	set(gca, 'Position', [0.045 0.015 0.36 0.48], 'FontSize', 12)
	subimage(uint8(pforecast.shadow{1}+1), naturalskymap)
	hold on
	scatter(site_X, site_Y, 'sk', 'SizeData', 100, 'LineWidth', 2)
	if nargin == 3
		text(site_X+25, site_Y-45, sprintf('%s', site_name), 'Color', 'k', 'FontWeight', 'bold')
	end
	hold off
	xlabel('Distance [m]')
	ylabel('Distance [m]')
	title('Shadow map over UCSD')
	
	% Zoomed in nowcast kt over DEMROES kt (bottom right subplot)
	subplot(2,3,[5 6], 'replace');
	set(gca, 'Position', [0.45 0.05 0.53 0.41], 'FontSize', 12)
	plot(timeUTC, nowcast_kt, 'r', 'LineWidth', 2);
	hold on
	plot(timeUTC, GHI_kt, 'g', 'LineWidth', 2);
	scatter(timeUTC(pidx), nowcast_kt(pidx), 'or', 'LineWidth', 2)
	hold off
	axis tight
	xlim([timeUTC(pidx)-15/24/60 timeUTC(pidx)+15/24/60])
	datetick('x', 'keeplimits')
	box on
	xlabel('Time (UTC) [hh:mm]')
	ylabel('Clear Sky Index kt [-]')
	legend('Nowcast kt', sprintf('%s DEMROES kt', site_name), 'Location', 'SouthWest')
	ktplotx = get(gca, 'xlim');
	ktploty = get(gca, 'ylim');
	if isfield(pcloudheight, 'metar')
		heightvar = pcloudheight.metar;
		heighttype = 'METAR';
	elseif isfield(pcloudheight, 'Ceilometer')
		heightvar = pcloudheight.Ceilometer;
		heighttype = 'CEILO';
	end
	velheading = acosd(dot([pcloudmotion.u.sky(pidx) pcloudmotion.v.sky(pidx)],[0 1])/norm([pcloudmotion.u.sky(pidx) pcloudmotion.v.sky(pidx)]));
	text(ktplotx(2) - 8.2/24/60, ktploty(1) + (ktploty(2) - ktploty(1))*0.08, sprintf('%s Cloud Height: %.0f m\nVelocity: %.1f m/s @ %.0f%c', heighttype, heightvar, norm([pcloudmotion.u.sky(pidx) pcloudmotion.v.sky(pidx)]), velheading, char(176)), 'FontSize', 12)
	title(sprintf('Nowcast kt overlaid on DEMROES calculated kt for %s (UTC)', datestr(timeUTC(pidx), 'mmm-dd-yyyy HH:MM:SS')))
		
	% Write everything to video
	currFrame=getframe(b);
	writeVideo(vidObj,currFrame);
end
	
% Close video Object
close(vidObj)

end