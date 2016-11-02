function smap(forecast)

% DEMROES locations
DEMROES.X = [647 679 764 2 1155 691];
DEMROES.Y = [759 545 514 1041 641 539];

%%% Copy & Paste from getGHI.m %%%

outputDir = sprintf('C:/Users/Handa/Documents/_My_Files/Videos/Smaps/%s', datestr(forecast.time(1),'yyyymmdd'));
if(exist(outputDir,'dir') == 0)
	mkdir(outputDir);
end

javaCheckInit();
target = siDeployment('UCSD');

for i = 1:31
	timepst(i) = addtodate(forecast.time(i), -8, 'hour');
end

% Get day of year
StartDate = datenum(forecast.time(1));
PrevYear = datenum(str2double(datestr(StartDate, 'yyyy'))-1, 12,31);
DOY = floor(StartDate-PrevYear);

% Load GHI data
for ID = 1:length(target.footprint.GHInames)
	GHIdata(ID) = load(siNormalizePath(sprintf('%s/%s/%s_%i_%i.mat','$KleisslLab4TB1/database/DEMROES/1s_by_DoY', target.footprint.GHInames{ID}, target.footprint.GHInames{ID}, str2double(datestr(forecast.time(1), 'yyyy')), DOY)));
	GHIstart = find(abs(GHIdata(ID).time_day(:,1) - addtodate(forecast.time(1), -8, 'hour')) < datenum([0000 00 00 00 00 01])); % Search for start time in DEMROES data and convert into UTC time. DEMROES -8 hours w.r.t. UTC.
	% Concatenate
	GHI(:,ID) = GHIdata(ID).GHI_day(GHIstart:GHIstart+900);
end

time = zeros(length(forecast.time),1);
time(1) = addtodate(forecast.time(1), -8, 'hour');
timeUTC(1) = forecast.time(1);

for i = 2:length(GHI(:,1))
	time(i) = addtodate(time(i-1), 1, 'second');
	timeUTC(i) = addtodate(timeUTC(i-1), 1, 'second');
end

GHIavg = zeros(length(time),1);
idx = 0;
for n = 1:length(GHI(:,1))
	idx = idx+1;
	GHIavg(idx) = mean(GHI(n,:));
end

% Convert DEMROES GHI Data to kt

csk = clearSkyIrradiance( target.ground.position , timeUTC, target.tilt, target.azimuth );

GHI_kt = GHIavg./csk.gi;

%%% End Copy & Paste

name = datestr(forecast.time(1));
vidObj = VideoWriter(sprintf('%s/%s-%s-%s-%s.avi',outputDir, name(1:11),name(13:14),name(16:17),name(19:20)));
vidObj.Quality = 100;
vidObj.FrameRate=5;
open(vidObj);

smapcm = [ 0 0 0; 1 1 1; 1 1 0; 0 0 1];

a = figure;
set(a, 'Position', [200 200 1120 420])
ticks = timepst(1):1/24/60:timepst(31);
ticklabels = cell(size(ticks));
for t = 1:length(ticks)
	ticklabels{t} = ticks(t);
end
h(1:2) = [0 0];
for i = 1:31
	% Shadow map
	sub1 = subplot(1,2,1);
	set(sub1, 'Position', [0.13 0.11 0.37 0.815])
	date = datestr(forecast.time(i));
	imagesc(forecast.shadow{i}, [-1 2])
	hold on
	scatter(DEMROES.X, DEMROES.Y, 'x', 'r')
	ylabel('Distance [m]')
	xlabel('Distance [m]')
	title(sprintf('Shadow map for %s', datestr(forecast.time(i))))
	colorbar
	colormap(smapcm)
	hold off
	text(1300, 1000,sprintf('%s:%s:%s', date(13:14), date(16:17), date(19:20)))
		
	% GHI
	subplot(1,2,2)
	h(2) = plot(timepst, forecast.ktavg,'r');
	hold on
	scatter(timepst(i), forecast.ktavg(i), 'm')
	set(gca, 'XTick', ticks, 'XTickLabel', ticklabels)
	datetick('x')
	xlim([timepst(1) timepst(31)])
	h(1) = plot(time, GHI_kt);
	xlabel('Time (PST) [hh:mm]')
	ylabel('Clear Sky Index (kt)')
	title(sprintf('GHI Forecast & DEMROES data for %s to %s (PST)', datestr(timepst(1)), datestr(timepst(length(timepst)))))
	legend(h, 'Avg DEMROES calculated kt', 'Forecast kt', 'Location', 'Best')
	hold off
	currFrame=getframe(a);
	writeVideo(vidObj,currFrame);
	
end

close(vidObj)

end