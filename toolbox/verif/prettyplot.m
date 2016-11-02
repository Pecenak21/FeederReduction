function [out] = prettyplot(imager, timestart, timeend, horizon, plotname, site_ID, option)
% prettyplot(imager, [timestart], [timeend], [horizon], [plotname], [siteID], [options...]);
% This function will create different plots based on forecast output data
% and real time reading from irradiance sensors or solar panel output.
% This function is currently tailored to work with DEMROES data at UCSD. If
% another site is set up, the code needs to be modified accordingly (Check on
% TODO items in this file to find out what needs to be changed).
%
% Input:
%			imager:		siImager or imager object for sky imager
%			timestart:  (optional) datenum or time string (e.g. 'yyyymmddHHMMSS') for start time (in UTC). Default: timeend - 30 minutes
%			timeend:	(optional) datenum or time string for end time (in UTC). Default: nowUTC
%			horizon:	(optional) array of forecast horizons in minue for plotting the forecast plot. Default: [0 5 10 15]
%			plotname:	(optional) name for plot of interest. Default: 'default'. Look below for more options.
%			site_ID:    (optional) ID for site of interest. NOTE: only works for DEMROES data at UCSD for the moment. Look below for more details.
%			option:     (optional) more options including 'fontsize','dirpath', 'quiet', 'width', 'height', 'visible','title','cache'
%a
% Examples of use:
%			i2 = siImager('USI_1_2'); % create imager
%			stime = datenum([2012 12 07 19 00 00]);
%			etime = datenum([2012 12 07 19 30 00]);
%			prettyplot(i2, stime, etime, [0 1 5], 'default', 3);	% Uses DEMROES data from EBU2 and forecast data using only EBU2's footprint.
%			prettyplot(i2);										% plot for last 30 minutes
%			prettyplot(i2, stime, etime, [], 'shadowmap');		% plot shadowmap with default time horizon
%
% Plot Names:
% 'default' - Plots the 'default' plot with the image, cloud decision, advection and forecast in the same plot
% 'shadowmap' - Only plots shadow map
% 'forecast' - Only plots forecast time series overlaid on DEMROES data
% 'video' - Makes a video of the 'default' plot with the forecast time series
%			graph zoomed in, showing data +- 15 minutes from the time of
%			the current frame. Also includes cloud height, cloud speed,
%			and heading as a text readout.
%
% Site ID's (only for UCSD):
% 1: BMSB, 2: CMRR, 3: EBU2, 4: HUBB, 5: MOCC, 6: POSL
%
% If no site ID is specified, prettyplot will average DEMROES data across
% all DEMROES stations, and the forecasts will be for the entire UCSD
% footprint.
%
% IMPORTANT NOTE REGARDING DATA LOCATION: if you want to set your forecast
% path differently from the one in forecast.conf, set option.dirpath to
% your interested output folder, in which each imager has its own folder
% named after its name. Default is the outputDir param from forecast.conf.

%% Process important inputs
if ~exist('imager','var') || ~(isa(imager,'imager') || ischar(imager))
	error('prettyplot:invalidInput','Please use a valid ''Imager'' object as input! Duhhh!');
end
if(ischar(imager))
	imager = siImager(imager);
end
if ~exist('timeend','var') || isempty(timeend)
	timeend = nowUTC; % default to now
else
	if ischar(timeend)
		timeend = datenum(timeend, 'yyyymmddHHMMSS');
	end
end
if ~exist('timestart','var') || isempty(timestart)
	timestart = timeend - 30/60/24; % default to 30 minutes from timeend
else
	if ischar(timestart)
		timestart = datenum(timestart, 'yyyymmddHHMMSS');
	end
end
if ~exist('horizon','var') || isempty(horizon)
	horizon = [0 5 10 15]; % in minutes
end
if ~exist('plotname','var') || isempty(plotname)
	plotname = 'default';
end

%% Caching to speed up performance
persistent target cache fcdata GHI_kt timeUTC timeinterval site_name site_X site_Y;

%% General parameter setup
% max forecast time ahead
max_time_forecast = 15; % in minutes

% interval between 2 forecasts
tinterval = 30; % in seconds

% Load config
conf = readConf(siGetConfPath('forecast.conf'));
if isinf(conf.endTime)
	ignoreSeriesFile = 1;
else
	ignoreSeriesFile = 0;
end

% Time start to time end in 30 second intervals
timeinterval = timestart:tinterval/24/3600:timeend;

%% handle optional inputs
%default value
opn = {'fontsize','quiet','dirpath','width','height','visible','title','cache','output','embed', 'save'};
op.fontsize = 12;
op.quiet = 0;
op.dirpath = conf.outputDir;
op.width = 1500;
op.height = 900;
op.visible = 1;
op.title = 1;
op.cache = 0;
op.output = 0;
op.embed = 0;
op.save = 0;

if ~exist('option','var') || isempty(option) || ~isstruct(option)
	option = op;
else
	for fi = opn
		if ~isfield(option,fi{1})
			option.(fi{1}) = op.(fi{1});
		end
	end
	id = ~ismember(fieldnames(option),opn);
	if sum(id) > 0
		fn = fieldnames(option);
		s = '';
		for i = find(id)
			s = sprintf('%s\n\t%s',s,fn{i});
		end
		warning('prettyplot:unsupportedoptions','There are %d unsupported options. They are:%s',sum(id),s);
	end
end

%% Just plot and skip the calculation if the setup has not changed
if op.cache
	if ~isempty(cache) && strcmpi(imager.name,cache.name) && (timestart == cache.timestart) && (timeend == cache.timeend) && isequal(horizon, cache.horizon)
		if (~exist('site_ID','var') && ~isfield(cache,'site_ID')) || (exist('site_ID','var') && isfield(cache,'site_ID') && strcmp(site_ID,cache.site_ID) || (isempty(site_ID) && isempty(cache.site_ID)) )
			plotnow(0);
			return;
		end
	end
end

%% Find interested directories
% Assign directories

%% Setup parameters specific to deployment site
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START SPECIFIC SITE SETUP %%%%%%%%%%%%%%%%%%
% TODO: will need to modify these parameters if a different site is used
switch(lower(conf.deployment))
	case 'ucsd'
		[site(1:6).name] = deal('BMSB','CMRR','EBU2','HUBB','MOCC','POSL');
		if isempty(target)
			target = siDeployment(conf.deployment);
			% Time Zone setup
		end
		if ~isfield(conf,'timezoneloc'), conf.timezoneloc = 'Los Angeles'; end
		if ~isfield(conf,'timezoneDST'), conf.timezoneDST = 0; end
		
		% If no site_ID specified, use average DEMROES data
		if ~exist('site_ID','var') || isempty(site_ID)
			% Average data
			site_start = 1;
			site_end = 6;
			site_name = 'Average';
			site_X = target.DEMROES_X;
			site_Y = target.DEMROES_Y;
			if ~option.quiet, fprintf('===%s Data Mode===\n', site_name); end
		else
			% Site specific data
			site_start = site_ID;
			site_end = site_ID;
			site_name = target.footprint.GHInames{site_start};
			site_X = target.DEMROES_X(site_start);
			site_Y = target.DEMROES_Y(site_start);
			if ~option.quiet, fprintf('===Site-specific Data Mode for %s===\n', site_name); end
		end
		
		%% Load DEMROES GHI data
		% change time to DST format (DEMROES format) to find the days and get correct data (data files are named based on DST format aaaaaaa)
		PrevYear = datevec(timeend);
		PrevYear = datenum(PrevYear(1)-1, 12,31);
		DOY = floor(toLocalTime(timestart,conf.timezoneloc,conf.timezoneDST)-PrevYear):...
			floor(toLocalTime(timeend,conf.timezoneloc,conf.timezoneDST)-PrevYear);
		
		% Initialize
		t = (timestart - 15/24/60 : 1/3600/24 : timeend)'; % time stamp with increment of 1 second; Read in GHI data 15 minutes prior (for error plot)
		GHI = nan( length(t), length(target.footprint.GHInames) );
		for id = site_start:site_end
			for d = DOY
				fn = siNormalizePath(sprintf('%s/%s/%s_%i_%i.mat','$KleisslLab4TB1/database/DEMROES/1s_by_DoY', target.footprint.GHInames{id}, target.footprint.GHInames{id}, str2double(datestr(timestart, 'yyyy')), d));
				if exist(fn,'file')
					% Load GHI data
					GHIdata = load(fn);
					GHIdata.time_day = toUTC( GHIdata.time_day , conf.timezoneloc, conf.timezoneDST );
					
					% Filter MOCC data by SZA only before Jan 18
					if( (id == 5) && (timestart < datenum([2013 01 18 08 00 00])) )
						MOCC_position.longitude = -117.222547;
						MOCC_position.latitude = 32.878443;
						MOCC_position.altitude = 23;
						MOCC_SUN = siSunPosition( GHIdata.time_day , MOCC_position );
						MOCC_maxzenith = atand(1.41 ./ cosd( abs(180 - [MOCC_SUN.azimuth]) ) );
							
						% only filter data before Jan 18
						GHIdata.GHI_day( (MOCC_SUN.zenith(:) > MOCC_maxzenith) & (GHIdata.time_day < datenum([2013 01 18 08 00 00])) ) = NaN;
					end
					% End MOCC filter
					
					% find matched time
					[lia,lib] = ismember( round(GHIdata.time_day.*(24*3600)), round( t.*(24*3600)) );
					if sum(lia) > 0
						lib(lib==0) = [];
						% assign corresponding GHI
						GHI(lib,id) = GHIdata.GHI_day( lia );
					end
				end
			end
		end
		
		% clean up GHI data by removing all timestamps that does not have any data from any station
		id = ~isnan(GHI);
		id = sum(id,2) > 0;
		% Average GHI data if needed
		if sum(id) > 0
			timeUTC = t(id);
			GHI = GHI(id,:);
			% refine one more time based on clear sky irradiance
			csk = clearSkyIrradiance( target.ground.position , timeUTC, target.tilt, target.azimuth );
			id = csk.gi > 0;
			timeUTC = timeUTC(id);
			GHI = GHI(id,:);
			% take average, change nan values to 0
			d = ~isnan(GHI);
			GHI(~d) = 0;
			GHIavg = sum(GHI,2)./sum(d,2);
			% Convert DEMROES GHI Data to kt
			GHI_kt = GHIavg./csk.gi(id);
			GHIclrsky = csk.gi(id);
			
			% apply a bit of quality control here
% 			GHI_kt(GHI_kt>1.5) = 1.5;
			
			GHI_kt_persistence = GHI_kt;
			GHIavg_persistence = GHIavg;
			timeUTC_persistence = timeUTC;
			% Delete beginning entries from timeUTC, GHIavg, GHI_kt so that
			% everything lines up with the start time
			idf = find(timeUTC == timestart);
			timeUTC(1:idf-1) = [];
			GHIavg(1:idf-1) = [];
			GHI_kt(1:idf-1) = [];
			
			% averaging
			for m = 1:floor(numel(GHIavg)/30)
				GHI30savg(m) = mean(GHIavg((m-1)*30+1:m*30));
				timeUTC30savg(m) = timeUTC(m*30);
			end
		else
			warning('prettyplot:nositedata','There is no measured data from this site: %s. Will skip plotting those lines on forecast plot.',conf.deployment);
		end
	case 'redlands'
		% TOBE FILLED IN
	case 'henderson'
		% TOBE FILLED IN
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END SPECIFIC SITE SETUP %%%%%%%%%%%%%%%%%%%%%%

%% Load forecast output files when needed for some specific plot
if( any( strcmp(plotname, {'forecast', 'forecastghi', 'forecastkt', 'default', 'rmse', 'rmseghi', 'video', 'ramprate', 'kt'}) ) )
	if max(horizon) < max_time_forecast
		max_time_forecast = max(horizon);
	end
	
	if ~exist('site_ID','var') || isempty(site_ID)
		data = loadforecastdata(imager,timestart-max_time_forecast/60/24,timeend,'forecast',{'time','ktavg','gi','kt','station'},option.dirpath, ignoreSeriesFile);
	else
		data = loadforecastdata(imager,timestart-max_time_forecast/60/24,timeend,'forecast',{'time','station'},option.dirpath, ignoreSeriesFile);
	end
	
	if isempty(data)
		warning('loadforecastdata:nodata','No data is available in specified time window!');
	else
		if ~option.quiet
			fprintf('Start time: %s UTC\n', datestr(timestart))
			fprintf('End time: %s UTC\n', datestr(timeend))
		end
		% Load data for corresponding forecast horizons
		fidx = 0;
		fcdata = struct();
		for fh = horizon
			fidx = fidx+1; % Counter
			% Label forecast horizon in fcdata structure (legend later)
			if fh == 0
				fcdata(fidx).forecastHorizon = 'USI nowcast (0 min)';
			else
				fcdata(fidx).forecastHorizon = sprintf('USI forecast (%i min)', fh);
			end
			
			offset = (max(max_time_forecast)*2 + 1) - fh*2;
			id = offset : size(data.time,1) - fh*2 ;
			fcdata(fidx).time = data.time(id,fh*2+1);
			if ~exist('site_ID','var') || isempty(site_ID)
				fcdata(fidx).kt = data.ktavg(id,fh*2+1);
				if isfield(data,'gi'), fcdata(fidx).gi = data.gi(id,fh*2+1); end
				if isfield(data,'station')
					for lid = 1:numel(id)
						for sid = 1:6
							fcdata(fidx).stationktavg(lid,sid) = data.station(id(lid),sid).ktavg(fh*2+1);
						end
					end
				end
			else
				fcdata(fidx).kt = data.station(id,fh*2+1);
			end
		end
		if ~isempty(fcdata)
			if length(fcdata) > 1
				forecast_time = fcdata(1).time;
			else
				forecast_time = fcdata.time;
			end
			if ~isempty(forecast_time), forecast_time = forecast_time(end); end
		else
			warning('prettyplot:noforecastdata','no forecast data is available in specified time');
		end
	end
end
%% caching results
cache.name = imager.name;
cache.timestart = timestart;
cache.timeend = timeend;
cache.horizon = horizon;
if exist('site_ID','var')
	cache.site_ID = site_ID;
else
	if isfield(cache,'site_ID'), cache = rmfield(cache,'site_ID');end
end

%% Plot
if ~option.output, plotnow(1); else out = plotnow(1); end
if ~option.quiet, fprintf('Done.\n'); end

%% Internal function to generate Pretty Plot
% TODO: this plot function might also be site specific now. Will change
% when needed later
	function [out] = plotnow(reloadFlag)
		% Input:
		%			reloadFlag		: trigger to reload the projection and forecast files when setup is changed
		
		if ~exist('reloadFlag','var')
			reloadFlag = 1;
		end
		
		% Initialize figure window
		if option.visible
			b = figure('Color', [1 1 1]);
		else
			b = figure('Color', [1 1 1], 'visible', 'off');
		end
		set(b, 'Position', [50 50 option.width option.height]);
		set(0,'defaulttextfontsize',option.fontsize);
		set(0,'defaultaxesfontsize',option.fontsize);
		
		% Create & assign color map
		naturalskymap = [ 1 1 1; 0 .4 .8; .9 .9 .9; .7 .7 .7 ];
		
		if ~exist('forecast_time','var') || isempty(forecast_time)
			forecast_time = timeend;
		end
		
		persistent pprojection pforecast pcloudheight;
		% Load pertinent files for end time
		if isempty(pprojection) || reloadFlag
			[x, p] = siForecastIOExists([option.dirpath '/' imager.name '/' datestr(forecast_time,'yyyymmdd')],'projection',forecast_time);
			if x,pprojection=load(p, 'rbr', 'dec');end
		end
		if isempty(pforecast) || reloadFlag
			[x, p] = siForecastIOExists([option.dirpath '/' imager.name '/' datestr(forecast_time,'yyyymmdd')],'forecast',forecast_time);
			if x, pforecast = load(p, 'shadow'); end
		end
		if isempty(pcloudheight) || reloadFlag
			[x, p] = siForecastIOExists([option.dirpath '/' imager.name '/' datestr(forecast_time,'yyyymmdd')],'cloudheight',forecast_time);
			if x, pcloudheight = load(p);end
		end
		
		switch(lower(plotname))
			case 'default'
				%% -- Default Plot
				% Plot Raw, RBR, CD, Shadowmap images for end time. Plot Forecast Time
				% Series and RMSE for compiled fcdata.
				
				% Load RBR and cloud decision from projection file
				% IMPORTANT: If projection file structure changes in the future, this
				% will need to be altered to accommodate
				
				% Raw image (top left)
				subplot(2,3,1, 'replace');
				set(gca, 'Position', [0.015 0.5 0.3 0.47], 'FontSize', option.fontsize)
				img = imread([imager.imageDir(forecast_time) '/' datestr(forecast_time, 'yyyymmdd') '/' datestr(forecast_time, 'yyyymmddHHMMSS') '_prev.jpg']);
				imagesc(img)
				axis image
				axis off
				if option.title, title('Raw Image'); end
				
				% RBR image (top middle)
				subplot(2,3,2, 'replace');
				set(gca, 'Position', [0.34 0.5 0.33 0.47], 'FontSize', option.fontsize)
				h = imagesc(pprojection.rbr, [0.5 1.5]);
				colorbar
				axis image
				axis off
				if option.title, title('Red-Blue-Ratio'); end
				
				% Cloud decision (top right)
				subplot(2,3,3, 'replace');
				set(gca, 'Position', [0.675 0.5 0.33 0.47], 'FontSize', option.fontsize)
				subimage(uint8(pprojection.dec+1), naturalskymap)
				axis image
				axis off
				if option.title, title('Cloud Decision'); end
				
				c = get(h,'CData');
				mask = isnan(pprojection.rbr);
				c(repmat(mask,[1 1 3])) = 255;
				set(h,'CData',c);
				
				
				% Shadow map (bottom left)
				subplot(2,3,4)
				set(gca, 'Position', [0.045 0.015 0.36 0.48], 'FontSize', option.fontsize)
				subimage(uint8(pforecast.shadow{1}+1), naturalskymap)
				hold on
				scatter(site_X, site_Y, 'sk', 'SizeData', 100, 'LineWidth', 2)
				if nargin == 5
					text(site_X+25, site_Y-45, sprintf('%s', site_name), 'Color', 'k', 'FontWeight', 'bold')
				end
				hold off
				xlabel('Distance [m]')
				ylabel('Distance [m]')
				% Height printout
				if isfield(pcloudheight, 'metar')
					heightvar = pcloudheight.metar;
					heighttype = 'METAR';
				elseif isfield(pcloudheight, 'Ceilometer')
					heightvar = pcloudheight.Ceilometer;
					heighttype = 'Ceilometer';
				end
				if option.title, title( sprintf('Shadow map over UCSD. %s cloud height: %.0f m.', heighttype, heightvar) ); end
				
				% Zoomed in nowcast kt over DEMROES kt (bottom right subplot)
				subplot(2,3,[5 6], 'replace');
				set(gca, 'Position', [0.45 0.05 0.53 0.41], 'FontSize', option.fontsize)
				% Plot DEMROES GHI data
				if exist('GHI_kt','var') && ~isempty(GHI_kt)
					plot(timeUTC, GHI_kt, 'g', 'LineWidth', 2);
				end
				hold on
				% Plot forecast time series
				clist = ['k' 'b' 'm' 'c' 'r'];
				for pidx = 1:numel(horizon)
					if ~isempty(fcdata(pidx).kt)
						plot(fcdata(pidx).time, fcdata(pidx).kt, clist(pidx), 'LineWidth', 2);
					else
						warning('prettyplot:datanotavailable','Forecast data for %.0f minute ahead is not availble! Either check the data or change forecast horizon.',horizon(pidx));
					end
				end
				hold off
				axis tight
				xlim([timestart timeend])
				datetick('x', 'keeplimits')
				box on
				xlabel('Time (UTC) [hh:mm]')
				ylabel('Clear Sky Index kt [-]')
				legend(sprintf('%s DEMROES', site_name), fcdata.forecastHorizon, 'Location', 'NorthWest')
				title(sprintf('USI forecast kt overlaid on DEMROES calculated kt for %s (UTC)', datestr(timeend, 'mmm-dd-yyyy HH:MM:SS')));
				
			case 'raw'
				%% -- Raw Plot
				img = imread([imager.imageDir '/' datestr(forecast_time, 'yyyymmdd') '/' datestr(forecast_time, 'yyyymmddHHMMSS') '_prev.jpg']);
				imagesc(img)
				axis image
				axis off
				if option.title, title('Raw Image'); end
				
			case 'rbr'
				%% -- RBR Plot
				if ~isempty(pprojection)
					h = imagesc(pprojection.rbr, [.5 1.5]);
					colorbar
					axis image
					axis off
					if option.title, title('Red-Blue-Ratio'); end
					% need to mask out the background nicely. Do not know how
					% to do it now.
					% 				c = get(h,'CData');
					% 				mask = isnan(pprojection.rbr);
					% 				c(repmat(mask,[1 1 3])) = 255;
					% 				set(h,'CData',c);
				end
				
			case 'clouddecision'
				%% -- Cloud Decision Plot
				if ~isempty(pprojection)
					subimage(uint8(pprojection.dec+1), naturalskymap)
					axis image
					axis off
					if option.title, title('Cloud Decision'); end
				end
			
			case 'shadowmap'
				%% -- Shadow Map Plot
				if ~isempty(pforecast)
					subimage(uint8(pforecast.shadow{1}+1), naturalskymap)
					hold on
					scatter(site_X, site_Y, 'sk', 'SizeData', 100, 'LineWidth', 2)
					if nargin == 5
						text(site_X+25, site_Y-45, sprintf('%s', site_name), 'Color', 'k', 'FontWeight', 'bold')
					end
					option.type = 'cloud';
					gmapOverlay(pforecast.shadow{1}, [-117.2534 -117.2127], [32.8695 32.8928], .7, option);
					
					% Height printout
					if isfield(pcloudheight, 'metar')
						heightvar = pcloudheight.metar;
						heighttype = 'METAR';
					elseif isfield(pcloudheight, 'Ceilometer')
						heightvar = pcloudheight.Ceilometer;
						heighttype = 'Ceilometer';
					end
					if option.title, title( sprintf('Shadow map over UCSD. %s cloud height: %.0f m.', heighttype, pcloudheight.height) ); end
				end
				
			case 'shadowmapWithoutOverlay'
				%% -- Shadow Map Plot
				if ~isempty(pforecast)
					subimage(uint8(pforecast.shadow{1}+1), naturalskymap)
					hold on
					scatter(site_X, site_Y, 'sk', 'SizeData', 100, 'LineWidth', 2)
					if nargin == 5
						text(site_X+25, site_Y-45, sprintf('%s', site_name), 'Color', 'k', 'FontWeight', 'bold')
					end
					hold off
					xlabel('Distance [m]')
					ylabel('Distance [m]')
					if option.title, title('Shadow map over UCSD'); end
					% Height printout
					if isfield(pcloudheight, 'metar')
						heightvar = pcloudheight.metar;
						heighttype = 'METAR';
					elseif isfield(pcloudheight, 'Ceilometer')
						heightvar = pcloudheight.Ceilometer;
						heighttype = 'Ceilometer';
					end
					if option.title, title( sprintf('Shadow map over UCSD. %s cloud height: %.0f m.', heighttype, pcloudheight.height) ); end
				end
				
			case 'forecastkt'
				%% -- Forecast Plot
				[~, stdidx] = ismember( round(fcdata(1).time.*(24*3600)), round(timeUTC.*(24*3600)));
				GHI_kt_subsample = GHI_kt(stdidx);

				if ~isempty(fcdata)
					if exist('GHI_kt','var') && ~isempty(GHI_kt)
						if ~option.embed
							plot(fcdata(1).time - 8/24, GHI_kt_subsample, 'Color', [0 0.9 0], 'LineWidth', 2);
						else
							plot(option.embed, fcdata(1).time, GHI_kt_subsample, 'g', 'LineWidth', 2);
						end
					end
% 					if exist('GHI_kt','var') && ~isempty(GHI_kt)
% 						if ~option.embed
% 							plot(timeUTC, GHI_kt, 'g', 'LineWidth', 2);
% 						else
% 							plot(option.embed, timeUTC, GHI_kt, 'g', 'LineWidth', 2);
% 						end
% 					end
					hold on
					% Plot forecast time series
					clist = ['k' 'b' 'm' 'c' 'r'];
					for pidx = 1:numel(horizon)
						if isfield(fcdata,'kt') && ~isempty(fcdata(pidx).kt)
							if ~option.embed
								plot(fcdata(pidx).time - 8/24, fcdata(pidx).kt, clist(pidx), 'LineWidth', 2);
							else
								plot(option.embed, fcdata(pidx).time, fcdata(pidx).kt, clist(pidx), 'LineWidth', 2);
							end
						else
							warning('prettyplot:datanotavailable','Forecast data for %.0f minute ahead is not availble! Either check the data or change forecast horizon.',horizon(pidx));
						end
					end
					hold off
					if ~option.embed
						axis tight
% 						xlim([timestart-8/24 timeend-8/24])
						ylim([0 1.4])
						datetick('x', 'keeplimits')
						box on
						xlabel('Time (PST) [hh:mm]')
						ylabel('Clear Sky Index kt [-]')
						legend('measured', fcdata.forecastHorizon, 'Location', 'SouthEast', 'Orientation', 'horizontal')
						legend('boxoff')
						if option.title, title(sprintf('USI forecast kt overlaid on DEMROES calculated kt for %s (UTC)', datestr(timeend, 'mmm-dd-yyyy HH:MM:SS'))); end
					end
				end
				if option.save
					options.Format = 'png';
					plotDir = [option.dirpath '/' imager.name '/Plots/ForecastKT/' datestr(timestart, 'yyyymmdd') '/'];
					if(exist(plotDir,'dir') == 0)
						mkdir(plotDir);
					end
					hgexport(b, [plotDir datestr(timestart, 'yyyymmdd') '_' sprintf('%i', horizon) '.png'], options);
				end
				
			case 'forecastghi'
				%% Forecast GHI plot
				if ~isempty(fcdata)

					% Subsample GHI data
					[~, stdidx] = ismember( round(fcdata(1).time.*(24*3600)), round(timeUTC.*(24*3600)));
					GHIavg_subsample = GHIavg(stdidx);

					%%% 30s SUBSAMPLE DEMROES GHI %%%
					% Plot measured GHI
					if exist('GHIavg_subsample','var') && ~isempty(GHIavg)
						if ~option.embed
							plot(fcdata(1).time, GHIavg_subsample, 'g', 'LineWidth', 2);
						end
					end
					%%% 1s DEMROES GHI %%%
% 					% Plot measured GHI
% 					if exist('GHIavg','var') && ~isempty(GHIavg)
% 						if ~option.embed
% 							plot(timeUTC, GHIavg, 'g', 'LineWidth', 2);
% 						else
% 							plot(option.embed, timeUTC, GHIavg, 'g', 'LineWidth', 2);
% 						end
% 					end
					%%%% 30s AVG DEMROES GHI %%%
% 					% Plot measured GHI
% 					if exist('GHI30savg','var') && ~isempty(GHIavg)
% 						if ~option.embed
% 							plot(timeUTC30savg, GHI30savg, 'g', 'LineWidth', 2);
% 						else
% 							plot(option.embed, timeUTC30savg, GHI30savg, 'g', 'LineWidth', 2);
% 						end
% 					end
					hold on
					% Plot forecast time series
					clist = {'k' 'b' 'm' 'c' 'r'};
					for pidx = 1:numel(horizon)
						if isfield(fcdata,'gi') && ~isempty(fcdata(pidx).gi)
							if ~option.embed
								plot(fcdata(pidx).time, fcdata(pidx).gi, clist{pidx}, 'LineWidth', 2);
							else
% % % % % 								fcdata(pidx).gi(isnan(fcdata(pidx).gi)) = 0;
% % % % % 								area(option.embed, fcdata(pidx).time, fcdata(pidx).gi, 'FaceColor', [0 .45 .9], 'EdgeColor', [0 .45 .9]);
								overundertest = fcdata(1).gi - GHIavg_subsample;
								under = (overundertest <= 0);
								over = (overundertest >= 0);
								gray_nodata = zeros(size(GHIavg_subsample));
								gray_nodata(fcdata(1).gi > 0) = GHIavg_subsample(fcdata(1).gi > 0);
								blue_fc = zeros(size(GHIavg_subsample));
								blue_fc(under) = fcdata(1).gi(under);
								blue_fc(~under) = GHIavg_subsample(~under);
								blue_fc(fcdata(1).gi == 0) = GHIavg_subsample(fcdata(1).gi == 0);
								red_fc = zeros(size(GHIavg_subsample));
								red_fc(over) = fcdata(1).gi(over);
								red_fc(~over) = GHIavg_subsample(~over);

								axes(option.embed)
								patch([1:size(GHIavg_subsample) size(GHIavg_subsample):-1:1], cat(1, GHIavg_subsample, flipud(gray_nodata)), 2, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none')
								patch([1:size(GHIavg_subsample) size(GHIavg_subsample):-1:1], cat(1, GHIavg_subsample, flipud(blue_fc)), 2, 'FaceColor', [0 0 1], 'EdgeColor', 'none')
								patch([1:size(GHIavg_subsample) size(GHIavg_subsample):-1:1], cat(1, GHIavg_subsample, flipud(red_fc)), 2, 'FaceColor', [1 0 0], 'EdgeColor', 'none')
								plot(option.embed, [length(GHIavg_subsample)/3 length(GHIavg_subsample)/3], [0 2000], '--k', 'LineWidth', 2)
								plot(option.embed, [2*length(GHIavg_subsample)/3 2*length(GHIavg_subsample)/3], [0 2000], '--k', 'LineWidth', 2)
							end
						else
							warning('prettyplot:datanotavailable','Forecast data for %.0f minute ahead is not available! Either check the data or change forecast horizon.',horizon(pidx));
						end
					end
					
% % % % % % % 					% Plot 30s subsampled DEMROES GHI if embedding (for calendarplot)
% % % % % % % 					if exist('GHIavg_subsample','var') && ~isempty(GHIavg)
% % % % % % % 						if option.embed
% % % % % % % 							plot(option.embed, fcdata(1).time, GHIavg_subsample, 'k', 'LineWidth', 2);
% % % % % % % 						end
% % % % % % % 					end
					
					hold off
					if ~option.embed
						axis tight
						xlim([timestart timeend])
						datetick('x', 'keeplimits')
						box on
						xlabel('Time (UTC) [hh:mm]')
						ylabel('Global Horizontal Irradiance [W m^{-2}]')
						legend(sprintf('%s DEMROES', site_name), fcdata.forecastHorizon, 'Location', 'NorthWest')
						legend('boxoff')
						if option.title, title(sprintf('USI forecast GHI overlaid on DEMROES calculated GHI for %s (UTC)', datestr(timeend, 'mmm-dd-yyyy HH:MM:SS'))); end
					end
					if option.save
						options.Format = 'png';
						plotDir = [option.dirpath '/' imager.name '/Plots/ForecastGHI/' datestr(timestart, 'yyyymmdd') '/'];
						if(exist(plotDir,'dir') == 0)
							mkdir(plotDir);
						end
						hgexport(b, [plotDir datestr(timestart, 'yyyymmdd') '_' sprintf('%i', horizon) '.png'], options);
					end
				end
				
			case 'forecast'
				if ~isempty(fcdata)
					%% will plot 3 horizontal subplots: the GHI time series, difference in GHIs of forecasts and real , and error.
					subplot(3,1,1);
					hold on
					if exist('GHIavg','var') && ~isempty(GHIavg)
						if ~option.embed
							plot(timeUTC, GHIavg, 'k', 'LineWidth', 1);
						else
							plot(option.embed, timeUTC, GHIavg, 'k', 'LineWidth', 2);
						end
					end
					clist = ['b' 'm' 'g' 'c' 'r'];
					for pidx = 1:numel(horizon)
						if isfield(fcdata,'gi') && ~isempty(fcdata(pidx).gi)
							if ~option.embed
								plot(fcdata(pidx).time, fcdata(pidx).gi, clist(pidx), 'LineWidth', 2);
							else
								plot(option.embed, fcdata(pidx).time, fcdata(pidx).gi, clist(pidx), 'LineWidth', 2);
							end
						else
							warning('prettyplot:datanotavailable','Forecast data for %.0f minute ahead is not availble! Either check the data or change forecast horizon.',horizon(pidx));
						end
					end

					if ~option.embed
						axis tight
						xlim([timestart timeend])
						ylim([0 1000]);
						set(gca,'YTick',0:200:1000)
						grid on
						box on
						ylabel('GHI [W/m^{2}]')
						datetick('x', 'keeplimits')
						set(gca,'xticklabel',[])
						legend(sprintf('%s DEMROES', site_name), fcdata.forecastHorizon, 'Location', 'NorthEast')
						if option.title, title(sprintf('USI forecast GHI overlaid on DEMROES calculated GHI for %s (UTC)', datestr(timeend, 'mmm-dd-yyyy HH:MM:SS'))); end
					end

					% subplot2: diff in GHI
					subplot(3,1,2);
					hold on
					clist = ['b' 'm' 'g' 'c' 'r'];
					if exist('GHIavg','var') && ~isempty(GHIavg)
						for pidx = 1:numel(horizon)
							if isfield(fcdata,'gi') && ~isempty(fcdata(pidx).gi)
								% find matching
								[y, x] = ismember(fcdata(pidx).time,timeUTC);
								if ~option.embed
									plot(fcdata(pidx).time(y), fcdata(pidx).gi(y)-GHIavg(x(x>0)), clist(pidx), 'LineWidth', 2);
								else
									plot(option.embed, fcdata(pidx).time(y), fcdata(pidx).gi(y)-GHIavg(x(x>0)), clist(pidx), 'LineWidth', 2);
								end
							else
								warning('prettyplot:datanotavailable','Forecast data for %.0f minute ahead is not availble! Either check the data or change forecast horizon.',horizon(pidx));
							end
						end
					end
					if ~option.embed
						axis tight
						xlim([timestart timeend])
						ylim([-500 500])
						grid on
						set(gca,'YTick',[-400:200:400]);
						box on
						ylabel('GHI_{fc}-GHI_{avg DEMROES} [W/m^{2}]')
						datetick('x', 'keeplimits')
						set(gca,'xticklabel',[])
						if option.title, title(sprintf('USI forecast GHI overlaid on DEMROES calculated GHI for %s (UTC)', datestr(timeend, 'mmm-dd-yyyy HH:MM:SS'))); end
					end

					% subplot3: error wrt to clear sky GHI
					subplot(3,1,3);
					hold on
					clist = ['b' 'm' 'g' 'c' 'r'];
					if exist('GHIavg','var') && ~isempty(GHIavg)
						for pidx = 1:numel(horizon)
							if isfield(fcdata,'gi') && ~isempty(fcdata(pidx).gi)
								% find matching
								[y, x] = ismember(fcdata(pidx).time,timeUTC);
								if ~option.embed
									plot(fcdata(pidx).time(y), abs((fcdata(pidx).gi(y)-GHIavg(x(x>0)))./GHIclrsky(x(x>0)))*100, clist(pidx), 'LineWidth', 2);
								else
									plot(option.embed, fcdata(pidx).time(y), abs((fcdata(pidx).gi(y)-GHIavg(x(x>0)))./GHIclrsky(x(x>0)))*100, clist(pidx), 'LineWidth', 2);
								end
							else
								warning('prettyplot:datanotavailable','Forecast data for %.0f minute ahead is not availble! Either check the data or change forecast horizon.',horizon(pidx));
							end
						end
					end
					if ~option.embed
						axis tight
						xlim([timestart timeend])
						ylim([0 100])
						grid on
						set(gca,'YTick',0:20:100)
						datetick('x', 'keeplimits')
						box on
						xlabel('Time (UTC) [hh:mm]')
						ylabel('|GHI_{fc}-GHI_{avg DEM}| / GHI_{clrSky} [%]')
						if option.title, title(sprintf('USI forecast GHI overlaid on DEMROES calculated GHI for %s (UTC)', datestr(timeend, 'mmm-dd-yyyy HH:MM:SS'))); end
					end
					samexaxis('abc','xmt','off','ytac','join','yld',1,'spaceRatio',.02)
				end
				
			case 'error'
				%%
				% This plot will show the accuracy of last hour and forecast for next 15 minutes
				if ~isempty(fcdata)
					subplot(2,1,1);
					if exist('GHI_kt','var') && ~isempty(GHI_kt)
						if ~option.embed
							plot(timeUTC, GHI_kt, 'g', 'LineWidth', 2);
						else
							plot(option.embed, timeUTC, GHI_kt, 'g', 'LineWidth', 2);
						end
					end
					hold on
					% Plot forecast time series
					clist = ['k' 'b' 'm' 'c' 'r'];
					for pidx = 1:numel(horizon)
						if isfield(fcdata,'kt') && ~isempty(fcdata(pidx).kt)
							if ~option.embed
								plot(fcdata(pidx).time, fcdata(pidx).kt, clist(pidx), 'LineWidth', 2);
							else
								plot(option.embed, fcdata(pidx).time, fcdata(pidx).kt, clist(pidx), 'LineWidth', 2);
							end
						else
							warning('prettyplot:datanotavailable','Forecast data for %.0f minute ahead is not availble! Either check the data or change forecast horizon.',horizon(pidx));
						end
					end
					hold off
					if ~option.embed
						axis tight
						xlim([timestart timeend])
						datetick('x', 'keeplimits')
						box on
						xlabel('Time (UTC) [hh:mm]')
						ylabel('Clear Sky Index kt [-]')
						legend(sprintf('%s DEMROES', site_name), fcdata.forecastHorizon, 'Location', 'NorthWest')
						if option.title, title(sprintf('USI forecast kt overlaid on DEMROES calculated kt for %s (UTC)', datestr(timeend, 'mmm-dd-yyyy HH:MM:SS'))); end
					end
				end
				
			case 'rmse'
				%% -- RMSE Plot
				% Calculate RMSE & MAE for each forecast horizon
				% Find idices of GHI data which match time of data
				% points and pull corresponding measured values to use
				% as standard by which forecast values are compared
				[~, stdidx] = ismember( round(fcdata(1).time.*(24*3600)), round(timeUTC.*(24*3600)));
				GHI_kt_std = GHI_kt(stdidx);
                beathorizon = NaN;
				for q = 1:length(horizon)
					% Compute persistence forecast
					[~, persistidx] = ismember( round((fcdata(q).time - horizon(q)/24/60).*(24*3600)), round(timeUTC_persistence.*(24*3600)));
					persistence_forecast = GHI_kt_persistence(persistidx);
					% Compute RMSE & MAE
					fcdata(q).RMSE = ( sqrt( nanmean( ( fcdata(q).kt(:) - GHI_kt_std(:) ).^2 ) ) / nanmean( GHI_kt_std(:) ) ) * 100;
					fcdata(q).pRMSE = ( sqrt( nanmean( ( persistence_forecast(:) - GHI_kt_std(:) ).^2 ) ) / nanmean( GHI_kt_std(:) ) ) * 100;
					fcdata(q).MAE = ( nanmean( abs( fcdata(q).kt(:) - GHI_kt_std(:) ) ) / nanmean( GHI_kt_std(:) ) ) * 100;
					fcdata(q).MBE = ( nanmean( fcdata(q).kt(:) - GHI_kt_std(:) ) / nanmean( GHI_kt_std(:) ) ) * 100;
					fcdata(q).pMAE = ( nanmean( abs( persistence_forecast(:) - GHI_kt_std(:) ) ) / nanmean( GHI_kt_std(:) ) ) * 100;
					fcdata(q).skill = 1 - ( fcdata(q).MAE/fcdata(q).pMAE );
					fcdata(q).RMSEskill = 1 - ( fcdata(q).RMSE/fcdata(q).pRMSE );
					if isnan(beathorizon) && ( fcdata(q).skill >= 0 )
						beathorizon = horizon(q);
					end
					
					% Compute average number of stations covered
					numstations = mean(~isnan(fcdata(q).stationktavg), 2).*6;
					fcdata(q).avgstationscovered = mean(numstations,1);
				end
				% Plot
				if ~option.output
					plot(horizon,cat(2,fcdata.RMSE), 'r', 'LineWidth', 3)
					hold on
					plot(horizon,cat(2,fcdata.MBE), 'b', 'LineWidth', 3)
					plot(horizon,cat(2,fcdata.MAE), 'Color', [0 .7 0], 'LineWidth', 3)
					plot(horizon,cat(2,fcdata.pMAE), 'LineStyle', '--', 'Color', [0 .7 0], 'LineWidth', 3)
					hold off
					title(sprintf('USI forecast (clear sky index kt) error metrics for %s to %s', datestr(timestart), datestr(timeend)))
					ylabel('Forecast Error [%]')
					xlabel('Forecast Horizon [minutes]')
					set(gca, 'XTick', horizon)
					legend('rRMSE', 'rMBE', 'rMAE', 'rMAE_{p}', 'Location', 'NorthEastOutside')
                    % Print forecast skills for 5, 10, 15 minutes
                    % Also print horizon at which we first beat persistence
                    x = xlim;
                    y = ylim;
                    text(1.01*x(2), 0.7*(y(2) - y(1)) + y(1), sprintf('5 minute skill: %.2f', fcdata(horizon == 5).skill))
                    text(1.01*x(2), 0.65*(y(2) - y(1)) + y(1), sprintf('10 minute skill: %.2f', fcdata(horizon == 10).skill))
                    text(1.01*x(2), 0.6*(y(2) - y(1)) + y(1), sprintf('15 minute skill: %.2f', fcdata(horizon == 15).skill))
                    if ~isnan(beathorizon)
						text(1.01*x(2), 0.55*(y(2) - y(1)) + y(1), sprintf('Skill > 0 at %d min', beathorizon ))
					end
				else
					% Dump forecast data into out variable
					out = fcdata;
				end
				
				if option.save
					options.Format = 'png';
					plotDir = [option.dirpath '/' imager.name '/Plots/kt_Error_wrt_FCHorizon_RMSE_skills/'];
					if(exist(plotDir,'dir') == 0)
						mkdir(plotDir);
					end
					hgexport(b, [plotDir datestr(timestart, 'yyyymmdd') '.png'], options);
					if(exist([plotDir '/Data/'],'dir') == 0)
						mkdir([plotDir '/Data/']);
					end
					save([plotDir '/Data/' datestr(timestart, 'yyyymmdd') '.mat'], 'fcdata')
				end
				
			case 'rmseghi'
				%% Calculate RMSE & MAE for each forecast horizon
				%% -- RMSE GHI plot
				% Find idices of GHI data which match time of data
				% points and pull corresponding measured values to use
				% as standard by which forecast values are compared
				[~, stdidx] = ismember( round(fcdata(1).time.*(24*3600)), round(timeUTC.*(24*3600)));
				GHIavg_std = GHIavg(stdidx);
                beathorizon = NaN;
				for q = 1:length(horizon)
					% Compute persistence forecast
					[~, persistidx] = ismember( round((fcdata(q).time - horizon(q)/24/60).*(24*3600)), round(timeUTC_persistence.*(24*3600)));
					persistence_forecast = GHIavg_persistence(persistidx);
					% Compute RMSE & MAE
					fcdata(q).RMSE = ( sqrt( nanmean( ( fcdata(q).gi(:) - GHIavg_std(:) ).^2 ) ) / nanmean( GHIavg_std(:) ) ) * 100;
					fcdata(q).MAE = ( nanmean( abs( fcdata(q).gi(:) - GHIavg_std(:) ) ) / nanmean( GHIavg_std(:) ) ) * 100;
					fcdata(q).MBE = ( nanmean( fcdata(q).gi(:) - GHIavg_std(:) ) / nanmean( GHIavg_std(:) ) ) * 100;
					fcdata(q).pMAE = ( nanmean( abs( persistence_forecast(:) - GHIavg_std(:) ) ) / nanmean( GHIavg_std(:) ) ) * 100;
					fcdata(q).skill = 1 - ( fcdata(q).MAE/fcdata(q).pMAE );
					if isnan(beathorizon) && ( fcdata(q).skill >= 0 )
						beathorizon = horizon(q);
					end
				end
				% Plot
				if ~option.output
					plot(horizon,cat(2,fcdata.RMSE), 'r', 'LineWidth', 3)
					hold on
					plot(horizon,cat(2,fcdata.MBE), 'b', 'LineWidth', 3)
					plot(horizon,cat(2,fcdata.MAE), 'Color', [0 .7 0], 'LineWidth', 3)
					plot(horizon,cat(2,fcdata.pMAE), 'LineStyle', '--', 'Color', [0 .7 0], 'LineWidth', 3)
					hold off
					title(sprintf('USI forecast error metrics for %s to %s', datestr(timestart), datestr(timeend)))
					ylabel('Forecast Error [%]')
					xlabel('Forecast Horizon [minutes]')
					set(gca, 'XTick', horizon)
					legend('rRMSE', 'rMBE', 'rMAE', 'rMAE_{p}', 'Location', 'NorthEastOutside')
                    % Print forecast skills for 5, 10, 15 minutes
                    % Also print horizon at which we first beat persistence
                    x = xlim;
                    y = ylim;
                    text(1.01*x(2), 0.7*(y(2) - y(1)) + y(1), sprintf('5 minute skill: %.2f', fcdata(horizon == 5).skill))
                    text(1.01*x(2), 0.65*(y(2) - y(1)) + y(1), sprintf('10 minute skill: %.2f', fcdata(horizon == 10).skill))
                    text(1.01*x(2), 0.6*(y(2) - y(1)) + y(1), sprintf('15 minute skill: %.2f', fcdata(horizon == 15).skill))
                    if ~isnan(beathorizon)
						text(1.01*x(2), 0.55*(y(2) - y(1)) + y(1), sprintf('Skill > 0 at %d min', beathorizon ))
					end
				else
					% Dump forecast data into out variable
					out = fcdata;
				end
				
				if option.save
					% Save .csv
					if numel(horizon) ~= 16
						error('Dump data for 0 to 15 min forecast horizons! Horizon = 0:15;')
					end
					csvDir = [option.dirpath '/' imager.name '/CSV/'];
					if(exist(csvDir,'dir') == 0)
						mkdir(csvDir);
					end
					fname = [csvDir datestr(timestart, 'yyyymmdd') '_' datestr(timestart, 'HHMMSS') '_' datestr(timeend, 'HHMMSS') '.csv'];
					fid = fopen(fname, 'w+');
					header = 'MATLAB Serial Date Number, Measured GHI, 0 min forecast [W m^-2], 1 min forecast [W m^-2], 2 min forecast [W m^-2], 3 min forecast [W m^-2], 4 min forecast [W m^-2], 5 min forecast [W m^-2], 6 min forecast [W m^-2], 7 min forecast [W m^-2], 8 min forecast [W m^-2], 9 min forecast [W m^-2], 10 min forecast [W m^-2], 11 min forecast [W m^-2], 12 min forecast [W m^-2], 13 min forecast [W m^-2], 14 min forecast [W m^-2], 15 min forecast [W m^-2]';
					fprintf(fid, '%s\n', header);
% 					fclose(fid);
					
					% Concatenate data for .csv file
					csvData(:,1) = fcdata(1).time;
					csvData(:,2) = GHIavg_std;
					for csvIdx = 1:16
						csvData(:,csvIdx+2) = fcdata(csvIdx).gi;
					end

					% Write to .csv!
% 					dlmwrite(fname, csvData(:,1), 'roffset', 1, 'precision', '%.34f', '-append');
					fprintf(fid, '%.34f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n', csvData');
					fclose(fid);
					
										
					% Save PNG
					options.Format = 'png';
					plotDir = [option.dirpath '/' imager.name '/Plots/Error_wrt_FCHorizon/'];
					if(exist(plotDir,'dir') == 0)
						mkdir(plotDir);
					end
					hgexport(b, [plotDir datestr(timestart, 'yyyymmdd') '.png'], options);
					if(exist([plotDir '/Data/'],'dir') == 0)
						mkdir([plotDir '/Data/']);
					end
					save([plotDir '/Data/' datestr(timestart, 'yyyymmdd') '.mat'], 'fcdata')
				end
				
			case 'ramprate'
				%% -- Ramprate Plot
				% set the index for a single plot
				fci = 1;
				%avgtimes = [0.5 1 5 15 60];
				avgtimes = [1  60];
				% ramprate plots
				if(length(fcdata)>1)
					warning('prettyplot:ramprate:multiHorizon','Ramp Rate plot only supports a single horizon');
				end
				% just do one of them
				[kt t] = interpavg(GHI_kt, timeUTC, 0.5, fcdata(fci).time);
				for avi = 1:length(avgtimes)
					rr(avi,1) = ramprates(kt,t,timestart,timeend,avgtimes(avi));
					rr(avi,2) = ramprates(fcdata(fci).kt,t,timestart,timeend,avgtimes(avi));
				end
				% make the plot for this guy
				%figure; % something like this was needed when we were looping fci
				plotstring = '';
				legendtext = {};
				for avi = 1:numel(avgtimes)
					for di = 1:2
						plotstring = sprintf('%s, rr(%d,%d).df_x,rr(%d,%d).pdf', plotstring, avi, di, avi, di);
					end
					legendtext = [legendtext {sprintf('measured, %d minute ramp', avgtimes(avi)), sprintf('forecast, %d minute ramp', avgtimes(avi))}];
				end
				h = eval(['semilogy(' plotstring(2:end) ')']);
				set(h(2),'color',get(h(1),'color'),'linestyle','--')
				set(h(4),'color',get(h(3),'color'),'linestyle','--')
				legend(legendtext);
				title(sprintf('Ramp Rate Distributions, %s', fcdata(fci).forecastHorizon));
				xlabel('kt Ramp Size')
				ylabel('Relative probability')

			case 'video'
				%% -- Video
                videoDir = [option.dirpath '/' imager.name '/Videos/' datestr(timestart, 'yyyymmdd') '/' sprintf('%i', horizon) '/'];
				if(exist(videoDir,'dir') == 0)
					mkdir(videoDir);
				end
				
                framecount = 0; % Used for PNG output
				for frameidx = 1:numel(timeinterval)
					% Plot Raw, RBR, CD, Shadowmap images for end time. Plot Forecast Time
					% Series and RMSE for compiled fcdata.
					
					% Load RBR and cloud decision from projection file
					% IMPORTANT: If projection file structure changes in the future, this
					% will need to be altered to accommodate
					
					% Load pertinent files for end time
					ppname = sprintf('%s/%s/projection/projection_%s.mat', [option.dirpath '/' imager.name], datestr(timeinterval(frameidx), 'yyyymmdd'), datestr(timeinterval(frameidx), 'yyyymmddHHMMSS'));
					pfname = sprintf('%s/%s/forecast/forecast_%s.mat', [option.dirpath '/' imager.name], datestr(timeinterval(frameidx), 'yyyymmdd'), datestr(timeinterval(frameidx), 'yyyymmddHHMMSS'));
					pchname = sprintf('%s/%s/cloudheight/cloudheight_%s.mat', [option.dirpath '/' imager.name], datestr(timeinterval(frameidx), 'yyyymmdd'), datestr(timeinterval(frameidx), 'yyyymmddHHMMSS'));
					pcmname = sprintf('%s/%s/cloudmotion/cloudmotion_%s.mat', [option.dirpath '/' imager.name], datestr(timeinterval(frameidx), 'yyyymmdd'), datestr(timeinterval(frameidx), 'yyyymmddHHMMSS'));
					
					prev_cmname = sprintf('%s/%s/cloudmotion/cloudmotion_%s.mat', [option.dirpath '/' imager.name], datestr(timeinterval(frameidx), 'yyyymmdd'), datestr((timeinterval(frameidx) - horizon(1)/24/60), 'yyyymmddHHMMSS'));
					prev_name = sprintf('%s/%s/forecast/forecast_%s.mat', [option.dirpath '/' imager.name], datestr(timeinterval(frameidx), 'yyyymmdd'), datestr((timeinterval(frameidx) - horizon(1)/24/60), 'yyyymmddHHMMSS'));
					if exist(ppname,'file') && exist(pfname,'file') && exist(pchname,'file') && exist(pcmname,'file') && exist(prev_name,'file') && exist(prev_cmname, 'file')
						pprojection = load(ppname, 'rbr', 'dec', 'mask');
						pforecast = load(pfname, 'shadow', 'gi', 'time');
						pcloudheight = load(pchname);
						pcloudmotion = load(pcmname, 'u', 'v');
						
						prev_cm = load(prev_cmname, 'u', 'v');
						prev_forecast = load(prev_name, 'shadow');
					else
						continue
					end
					
					% Raw image (top left)
					subplot(2,3,1, 'replace');
					set(gca, 'Position', [0.015 0.5 0.3 0.47], 'FontSize', option.fontsize)
					img = imread(sprintf('%s/%s/%s_prev.jpg', siNormalizePath(imager.imageDir(timeinterval(frameidx))),datestr(timeinterval(frameidx), 'yyyymmdd'), datestr(timeinterval(frameidx), 'yyyymmddHHMMSS')));
					imagesc(img)
					axis image
					axis off
					title('Raw Image')
					
					% RBR image (top middle)
					subplot(2,3,2, 'replace');
					set(gca, 'Position', [0.34 0.5 0.33 0.47], 'FontSize', option.fontsize)
					mask = isnan(pprojection.rbr);
					h = imagesc(pprojection.rbr, [0.5 1.5]);
					colorbar
					axis image
					axis off
					title('Red-Blue-Ratio')
					
					% Cloud decision (top right)
					subplot(2,3,3, 'replace');
					set(gca, 'Position', [0.675 0.5 0.33 0.47], 'FontSize', option.fontsize)
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
					set(gca, 'Position', [0.045 0.015 0.36 0.48], 'FontSize', option.fontsize)
					subimage(uint8(prev_forecast.shadow{horizon(1)*2+1}+1), naturalskymap)
					hold on
					scatter(site_X, site_Y, 'sk', 'SizeData', 100, 'LineWidth', 2)
					if nargin == 5
						text(site_X+25, site_Y-45, sprintf('%s', site_name), 'Color', 'k', 'FontWeight', 'bold')
					end
					% Show cm arrow
					quiver(720, 480, prev_cm.u.sky(end)*30/2.5, -prev_cm.v.sky(end)*30/2.5, 0, 'LineWidth', 2, 'Color', 'k')
					hold off
					set(gca, 'XTickLabel', num2str(str2num(get(gca, 'XTickLabel')).*2.5))
					set(gca, 'YTickLabel', num2str(str2num(get(gca, 'YTickLabel')).*2.5))
					xlabel('Distance [m]')
					ylabel('Distance [m]')
					title('Shadow map over UCSD')
					
					% Zoomed in nowcast kt over DEMROES kt (bottom right subplot)
					subplot(2,3,[5 6], 'replace');
					set(gca, 'Position', [0.45 0.05 0.53 0.41], 'FontSize', option.fontsize)
% % % 					% Plot DEMROES kt data
% % % 					if exist('GHI_kt','var') && ~isempty(GHI_kt)
% % % 						plot(timeUTC, GHI_kt, 'g', 'LineWidth', 2);
% % % 					end
					% Plot DEMROES kt data
					if exist('GHIavg','var') && ~isempty(GHIavg)
						plot(timeUTC, GHIavg, 'g', 'LineWidth', 2);
					end
					hold on
					% Plot forecast time series
					clist = ['k' 'b' 'y' 'c' 'r'];
					for pidx = 1:numel(horizon)
						plot(fcdata(pidx).time, fcdata(pidx).gi, clist(pidx), 'LineWidth', 2);
					end
					axis tight
					
					% Plot 15 min forecast generated at current time
					plot(pforecast.time, pforecast.gi, '--r', 'LineWidth', 2)
					
					% Plot vertical dotted line indicating time
					% corresponding with current frame
					plot([timeinterval(frameidx) timeinterval(frameidx)], [0 2000], ':k')
					plot([pforecast.time(horizon(1)*2 + 1) pforecast.time(horizon(1)*2 + 1)], [0 2000], ':k') % Vertical line corresponding to current forecast for horizon
					hold off
					% Plot will be zoomed in to show +- 15 minutes of data
					xlim([timeinterval(frameidx)-15/24/60 timeinterval(frameidx)+15/24/60])
					ktplotx = get(gca, 'xlim');
					ktploty = get(gca, 'ylim');
					ylim([0 ktploty(2)+0.1])
					ktploty = get(gca, 'ylim');
					datetick('x', 'keeplimits')
					box on
					xlabel('Time (UTC) [hh:mm]')
					ylabel('Global Horizontal Irradiance [W m^{-2}]')
					legend(sprintf('%s actual', site_name), fcdata.forecastHorizon, 'Current forecast', 'Location', 'SouthWest', 'Orientation', 'horizontal')
                    legend('boxoff')
					%% Height, velocity printout %%%
					if isfield(pcloudheight, 'metar')
						heightvar = pcloudheight.metar;
						heighttype = 'METAR';
					elseif isfield(pcloudheight, 'Ceilometer')
						heightvar = pcloudheight.Ceilometer;
						heighttype = 'CEILO';
					end
					velheading = sign(pcloudmotion.u.sky(end)).*acosd(dot([pcloudmotion.u.sky(end) pcloudmotion.v.sky(end)],[0 1])/norm([pcloudmotion.u.sky(end) pcloudmotion.v.sky(end)]));
					if velheading < 0
						velheading = velheading + 360;
					end
					% Text in bottom right corner
					% 						text(ktplotx(2) - 8.2/24/60, ktploty(2)*0.07, sprintf('%s Cloud Height: %.0f m\nVelocity: %.1f m/s @ %.0f%c', heighttype, heightvar, norm([pcloudmotion.u.sky(end) pcloudmotion.v.sky(end)]), velheading, char(176)), 'FontSize', option.fontsize)
					% Text in top left corner
					text(ktplotx(1) + 0.4/24/60, ktploty(2)*0.93, sprintf('%s Cloud Height: %.0f m\nVelocity: %.1f m/s @ %.0f%c', heighttype, heightvar, norm([pcloudmotion.u.sky(end) pcloudmotion.v.sky(end)]), velheading, char(176)), 'FontSize', option.fontsize)
					%% Wrap up
% 					title(sprintf('USI forecast kt overlaid on DEMROES calculated kt for %s (UTC)', datestr(timeinterval(frameidx), 'mmm-dd-yyyy HH:MM:SS')));
					title(sprintf('USI forecast GHI overlaid on DEMROES GHI for %s (UTC)', datestr(timeinterval(frameidx), 'mmm-dd-yyyy HH:MM:SS')));
                    % Write everything to PNG
                    options.Format = 'png';
                    hgexport(b, [videoDir '/' sprintf('%04d', framecount) '.png'], options);
                    framecount = framecount + 1;
                end
                % Write video using avconv (only works on Linux)
                if isunix
					mp4Dir = [videoDir '/../../MP4/'];
					if(exist(mp4Dir,'dir') == 0)
						mkdir(mp4Dir);
					end
					mp4filename = sprintf('%s/%s_%s.mp4', mp4Dir, datestr(timeinterval(1),'yyyymmdd'), sprintf('%i', horizon));
					if(exist(mp4filename,'file'))
						error('The automatically chosen output file name already exists.  Cancelling!');
					end
% 					error('dailymovie:missingX264','dailymovie requires x264 to be installed.  try ''sudo apt-get install x264''.\nYou will also need to do something like ''ln -s /usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.16 /usr/local/MATLAB/R2012b/sys/os/glnxa64/libstdc++.so.6''');
					unix(['x264 ' videoDir '/%04d.png --crf 18 -o ' mp4filename]);
                end
			case 'kt'
					if ~isempty(fcdata)
					hold on
					plot(data.time(:,1), data.kt(:,1), 'k', 'LineWidth', 2)
					plot(data.time(:,1), data.kt(:,2), 'r', 'LineWidth', 2)
					plot(data.time(:,1), data.kt(:,3), 'b', 'LineWidth', 2)
					hold off
					if ~option.embed
						axis tight
						datetick('x', 'keeplimits')
						box on
						xlabel('Time (UTC) [hh:mm]')
						ylabel('kt [-]')
						legend('Thick', 'Thin', 'Clear', 'Location', 'NorthWest')
						title(sprintf('kt pdf peak values for %s (UTC)', datestr(data.time(1,1), 'mmm-dd-yyyy')))
					end
				end
		end
		
		hold off;
		
	end

end
