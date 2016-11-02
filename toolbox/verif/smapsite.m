function smapsite(forecast, ID)

% DEMROES locations
DEMROES.X = [647 679 764 2 1155 691];
DEMROES.Y = [759 545 514 1041 641 539];

%%% Copy & Paste from getGHI.m %%%

outputDir = sprintf('C:/Users/Handa/Documents/_My_Files/Videos/Smaps/Site_%i/%s', ID, datestr(forecast.time(1),'yyyymmdd'));
if(exist(outputDir,'dir') == 0)
	mkdir(outputDir);
end

javaCheckInit();
target = siDeployment('UCSD');

for q = 1:31
	timepst(q) = addtodate(forecast.time(q), -8, 'hour');
end

% Get day of year
StartDate = datenum(forecast.time(1));
PrevYear = datenum(str2double(datestr(StartDate, 'yyyy'))-1, 12,31);
DOY = floor(StartDate-PrevYear);

% Load GHI data for SPECIFIC SITE
	GHIdata(ID) = load(siNormalizePath(sprintf('%s/%s/%s_%i_%i.mat','$KleisslLab4TB1/database/DEMROES/1s_by_DoY', target.footprint.GHInames{ID}, target.footprint.GHInames{ID}, str2double(datestr(forecast.time(1), 'yyyy')), DOY)));
	GHIstart = find(abs(GHIdata(ID).time_day(:,1) - addtodate(forecast.time(1), -8, 'hour')) < datenum([0000 00 00 00 00 01])); % Search for start time in DEMROES data and convert into UTC time. DEMROES -8 hours w.r.t. UTC.
	% Concatenate
	GHI = GHIdata(ID).GHI_day(GHIstart:GHIstart+300); % +300 for 5 min horizon, +900 for 15 min horizon


time = zeros(length(forecast.time),1);
time(1) = addtodate(forecast.time(1), -8, 'hour');

timeUTC(1) = forecast.time(1);

for q = 2:length(GHI)
	time(q) = addtodate(time(q-1), 1, 'second');
	timeUTC(q) = addtodate(timeUTC(q-1), 1, 'second');
end

% Convert DEMROES GHI Data to kt

csk = clearSkyIrradiance( target.ground.position , timeUTC, target.tilt, target.azimuth );

GHI_kt = GHI./csk.gi;
%%% End Copy & Paste

name = datestr(forecast.time(1));
vidObj = VideoWriter(sprintf('%s/%s-%s-%s-%s.avi',outputDir, name(1:11),name(13:14),name(16:17),name(19:20)));
vidObj.Quality = 100;
vidObj.FrameRate=5;
open(vidObj);

smapcm = [ 0 0 0; 1 1 1; 1 1 0; 0 0 1];

a = figure;
set(a, 'Position', [200 200 1120 420])
% Trim for purposes of plotting--match e.g. 1:11 to q below

for q = 1:31
	% Shadow map
	sub1 = subplot(1,2,1);
	set(sub1, 'Position', [0.13 0.11 0.37 0.815])
	date = datestr(forecast.time(q));
	imagesc(forecast.shadow{q}, [-1 2])
	hold on
	scatter(DEMROES.X(ID), DEMROES.Y(ID), 'x', 'r', 'SizeData', 400, 'LineWidth', 2)
	ylabel('Distance [m]')
	xlabel('Distance [m]')
	title(sprintf('Shadow map for %s', datestr(forecast.time(q))))
	colorbar
	colormap(smapcm)
	hold off
	text(1300, 1000,sprintf('%s:%s:%s', date(13:14), date(16:17), date(19:20)))
		
	% GHI
	subplot(1,2,2)
	plot(time, GHI_kt)
	hold on
	plot(timepst, forecast.station(ID).ktavg,'r')
	scatter(timepst(q), forecast.station(ID).ktavg(q), 'm')
	hold off
	datetick('x', 'keeplimits', 'keepticks')
	xlim([timepst(1) timepst(31)])
	xlabel('Time (PST) [hh:mm]')
	ylabel('Clear Sky Index (kt)')
	title(sprintf('GHI Forecast & %s data for %s to %s (PST)', target.footprint.GHInames{ID}, datestr(timepst(1)), datestr(timepst(length(timepst)))))
	legend(sprintf('%s calculated kt', target.footprint.GHInames{ID}), 'Forecast kt', 'Location', 'Best')
	currFrame=getframe(a);
	writeVideo(vidObj,currFrame);
	
end

close(vidObj)

end