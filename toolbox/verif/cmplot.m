function output = cmplot(date, inputDir, outputDir)
% Generates cloud motion plot for all forecast files within a certain
% date's folders

% inputDir requires to be st like xxx/USI_name/20121211

%% Process input
if ~exist('inputDir','var') || isempty(inputDir) 
	error('InputDir is required');
elseif ~exist(inputDir,'dir')
	error('InputDir doesn''t exist!');
end

if ~exist('outputDir','var') || isempty(outputDir)
	outputDir = 'tmp';
end

if ~exist(outputDir,'dir')
	mkdir(outputDir);
end

% Force date to be a string for consistency
if ~ischar(date)
	date = num2str(date);
end

conf.FontSize = 16;

%% Initial file loading
% Last cloud motion file
cm_files = dir(sprintf('%s/cloudmotion/', inputDir));
load(sprintf('%s/cloudmotion/%s', inputDir, cm_files(numel(cm_files)).name), 'cldfraction', 'time', 'u', 'v');
% Forecast files
forecast_files = dir(sprintf('%s/forecast/', inputDir));
% Cloud Height files
height_files = dir(sprintf('%s/cloudheight/', inputDir));
% % Cloud Decision files
% cd_files = dir(sprintf('%s/clouddecision/', inputDir));

%% Compute matching and cap errors
time_nowcast = nan(numel(cm_files)-3,1);
cldheight = nan(numel(cm_files)-3,1);

matching.totalpfalse = 0;
persist.totalpfalse = 0;

if (cldfraction(1) < 0.05) || (cldfraction(1) > 0.95)
	% Set cloud speeds to NaN
	u.sky(1) = NaN;
	v.sky(1) = NaN;
end

for q = 3:numel(cm_files)-1
	tic
	if (cldfraction(q-1) <= 0.05) || (cldfraction(q-1) >= 0.95) % (q-1) corresponds to cloud motion files, which exist for the first timestep (cldfraction, time, u, v)
		if cldfraction(q-1) <= 0.05
			fprintf('Processing %ith file of %i...Cloud fraction below 5%%.\n', q, numel(cm_files)-1)
		elseif cldfraction(q-1) >= 0.95
			fprintf('Processing %ith file of %i...Cloud fraction above 95%%.\n', q, numel(cm_files)-1)
		end
		nowcast = load(sprintf('%s/forecast/%s', inputDir, forecast_files(q+1).name), 'time');
		load(sprintf('%s/cloudheight/%s', inputDir, height_files(q+1).name));
		cldheight(q-2) = NaN;
		time_nowcast(q-2) = nowcast.time(1);
		
		% Set cloud speeds to NaN
		u.sky(q-1) = NaN;
		v.sky(q-1) = NaN;
		
		% Set errors to NaN
		err.matching(q-2) = NaN;
		err.persist(q-2) = NaN;
		err.cap(q-2) = NaN;
	else
		% Load files
		fprintf('Processing %ith file of %i...\n', q, numel(cm_files)-1)
		nowcast = load(sprintf('%s/forecast/%s', inputDir, forecast_files(q+1).name), 'advect', 'time');
		forecast = load(sprintf('%s/forecast/%s', inputDir, forecast_files(q).name), 'advect', 'time');
		load(sprintf('%s/cloudheight/%s', inputDir, height_files(q+1).name));
		cldheight(q-2) = metar;
		% 	load(sprintf('%s/clouddecision/%s', inputDir, cd_files(q+1).name), 'alpha');
		
		% Search for index corresponding to nowcast in forecast struct
		z = find(forecast.time == nowcast.time(1));
		
		% Matching err
		fmap = forecast.advect{z};
		fmap(fmap == 1) = 2; % Thin & thick cloud pixels are all = 2
		nmap = nowcast.advect{1};
		nmap(nmap == 1) = 2; % Thin & thick cloud pixels are all = 2
		matching.map = double(fmap - nmap);
		% Crop via NaNs
		matching.map(forecast.advect{z} == -1) = NaN;
		matching.map(nowcast.advect{1} == -1) = NaN;
		matching.mask = ~isnan(matching.map);
		matching.ptotal = sum(matching.mask(:));
		matching.errs = matching.map;
		matching.errs(matching.errs == 0) = NaN;
		matching.mask = ~isnan(matching.errs);
		matching.pfalse = sum(matching.mask(:));
		matching.totalpfalse = matching.totalpfalse + matching.pfalse;
		
		err.matching(q-2) = (matching.pfalse/matching.ptotal)*100;
		
		% Cap err
		pfmap = forecast.advect{1};
		pfmap(pfmap == 1) = 2; % Thin & thick cloud pixels are all = 2
		persist.map = double(pfmap - nmap);
		% Crop via NaNs
		persist.map(forecast.advect{z} == -1) = NaN; % Mask out same area for persistence as for matching so ptotal is the same in both error calculations
		persist.map(nowcast.advect{1} == -1) = NaN;
		persist.mask = ~isnan(persist.map);
		persist.ptotal = sum(persist.mask(:));
		persist.errs = persist.map;
		persist.errs(persist.errs == 0) = NaN;
		persist.mask = ~isnan(persist.errs);
		persist.pfalse = sum(persist.mask(:));
		persist.totalpfalse = persist.totalpfalse + persist.pfalse;
		
		err.persist(q-2) =(persist.pfalse/persist.ptotal)*100;
		err.cap(q-2) = (err.matching(q-2)/err.persist(q-2))*100;
		
		time_nowcast(q-2) = nowcast.time(1);
	end
	toc
end

%% Compute mean metrics

output.meanmatching = nanmean(err.matching);
output.stdmatching = nanstd(err.matching);
output.totalsumcap = (matching.totalpfalse / persist.totalpfalse) * 100;
output.sumcaptest = ( nansum(err.matching) / nansum(err.persist) ) * 100;
output.meancap = nanmean(err.cap);
output.stdcap = nanstd(err.cap);
output.avgcldfraction = nanmean(cldfraction);
% Compute cloud speed (cm.speed.mean is incorrect)
cloudspeed = sqrt( u.sky(1:end).^2 + v.sky(1:end).^2 );
%
output.avgcloudspeed = nanmean(cloudspeed);
output.avgcloudheight = nanmean(cldheight);
output.avgpersist = nanmean(err.persist);

%% Plot cloud speed, cloud height, cloud fraction, matching error, cap error
a = figure('Color', [1 1 1], 'Position', [50 50 1500 900], 'Visible', 'off');
ylabpos = time_nowcast(1) - 20/24/60; % First time minus 20 minutes

axes('Position', [0.07 0.79 0.9 0.16], 'FontSize', conf.FontSize);
hold on
plot(time,u.sky, 'b', 'LineWidth', 2), plot(time, v.sky, 'r', 'LineWidth', 2)
datetick('x', 15)
axis tight
xlim([time_nowcast(1) time_nowcast(numel(time_nowcast))])
set(gca, 'XTickLabel', [])
ylim([-60 60])
ylabel('Cloud Speed [m ^{s-1}]', 'Position', [ylabpos 0 1])
box on, grid on
set(gca, 'YTick', [-60 -40 -20 0 20 40 60]);
legend('E-W', 'N-S', 'Location', 'NorthEast', 'Orientation', 'horizontal')
legend('boxoff')
title( sprintf('Various time series for %s', datestr(time_nowcast(1), 'mmm dd, yyyy')) )
hold off

axes('Position', [0.07 0.61 0.9 0.16], 'FontSize', conf.FontSize);
plot(time_nowcast,cldheight, 'LineWidth', 2)
datetick('x', 15)
axis tight
xlim([time_nowcast(1) time_nowcast(numel(time_nowcast))])
set(gca, 'XTickLabel', [])
ylim([0 8000])
box on, grid on
set(gca, 'YTick', [0 2000 4000 6000 8000]);
ylabel('Cloud Height [m]', 'Position', [ylabpos 4000 1])

axes('Position', [0.07 0.43 0.9 0.16], 'FontSize', conf.FontSize);
plot(time,cldfraction, 'LineWidth', 2)
datetick('x', 15)
axis tight
xlim([time_nowcast(1) time_nowcast(numel(time_nowcast))])
set(gca, 'XTickLabel', [])
ylim([0 1])
box on, grid on
set(gca, 'YTick', [0 0.25 0.5 0.75 1]);
ylabel('Cloud Fraction [-]', 'Position', [ylabpos .5 1])

axes('Position', [0.07 0.25 0.9 0.16], 'FontSize', conf.FontSize);
plot(time_nowcast, err.matching, 'LineWidth', 2)
datetick('x', 15)
axis tight
xlim([time_nowcast(1) time_nowcast(numel(time_nowcast))])
set(gca, 'XTickLabel', [])
ylim([0 40])
box on, grid on
set(gca, 'YTick', [0 10 20 30 40]);
ylabel('Matching error [%]', 'Position', [ylabpos 20 1])

axes('Position', [0.07 0.07 0.9 0.16], 'FontSize', conf.FontSize);
hold on
plot(time_nowcast, err.cap, 'LineWidth', 2)
plot(time_nowcast, 100, ':k', 'LineWidth', 2)
hold off
datetick('x', 15)
axis tight
xlim([time_nowcast(1) time_nowcast(numel(time_nowcast))])
ylim([0 300])
ylabel('Cap error [%]', 'Position', [ylabpos 150 1])
box on, grid on
set(gca, 'YTick', [0 50 100 150 200 250 300]);
xlabel('Time [HH:mm] (UTC)')

%% Save outputs
% Save figure to PNG
options.Format = 'png';
hgexport(a, [outputDir '/' date '.png'], options);

% Save output
save([outputDir '/' date '.mat'], 'output')
end