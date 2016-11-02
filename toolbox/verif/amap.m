function amap(USI_ID,date)

if ischar(date) == 0
	date = sprintf('%i',date);
end

inputDir = sprintf('C:/Users/Handa/Desktop/ForecastOut/USI_1_%i/%s/forecast/',USI_ID,date);

% Load all forecast .mat files in directory
files = dir(inputDir);

% Define color map (Currently for -3 to 6)
ovcmap = [0 0 0; 0 0 0; 0 0 0; 1 1 1; 1 1 0; 1 0 0; 0 1 1; 0 0 0; 0 0 0; 0 0 1];

outputDir = './_My_Files/Videos/Amaps/';
if(exist(outputDir,'dir') == 0)
	mkdir(outputDir);
end

% Open video object
vidObj = VideoWriter(sprintf('%s%s-%s-%s.avi', outputDir, date(5:6), date(7:8), date(1:4)));
vidObj.Quality = 100;
vidObj.FrameRate = 3;
open(vidObj);

for i = 3:length(files)-1
	forecast = load(sprintf('%s/%s',inputDir,files(i).name));
	nowcast = load(sprintf('%s/%s',inputDir,files(i+1).name));

	figure
	% Plot nowcast as solid base
	imagesc(nowcast.advect{1}, [-3 6]);
	colormap(ovcmap)
	hold on
	
	% Search for index corresponding to nowcast in forecast struct
	q = find(forecast.time == nowcast.time(1));
	
	% Plot forecast as overlay -- set AlphaData to make everything
	% transparent except clouds
	forecastadvect = imagesc(forecast.advect{q}.*3, [-3 6]);
	ovalpha = get(forecastadvect, 'CData');
	ovalpha(ovalpha <= 0 | ovalpha == 4 | ovalpha == 5) = 0;
	ovalpha(ovalpha == 1 | ovalpha == 2 | ovalpha == 3 | ovalpha == 6) = 0.4;
	colormap(ovcmap)
	set(forecastadvect, 'AlphaData', ovalpha)
	axis image
	axis off
	title(sprintf('30s Forecast Cloud Map overlaid on Nowcast Cloud Map [%s (UTC)]', datestr(nowcast.time(1))))
	hold off
	set(gcf,'Position',[520 210 640 620])
	currFrame=getframe(gcf);
	writeVideo(vidObj,currFrame);
	close
end
close(vidObj)
end	