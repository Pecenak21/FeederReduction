function amap30to5(USI_ID,date,startstr)

if ischar(date) == 0
	date = sprintf('%i',date);
end
if ~ischar(startstr)
	startstr = num2str(startstr);
end

inputDir = sprintf('C:/Users/Handa/Desktop/ForecastOut/USI_1_%i/%s/forecast/',USI_ID,date);

% Load all forecast .mat files in directory
files = dir(inputDir);

for z = 1:length(files)
	match = strfind(files(z).name, startstr);
	if ~isempty(match)
		file_id = z; break
	end
end

% Define color map (Currently for -3 to 6)
ovcmap = [0 0 0; 0 0 0; 0 0 0; 1 1 1; 1 1 0; 1 0 0; 0 1 1; 0 0 0; 0 0 0; 0 0 1];

outputDir = './_My_Files/Videos/Amaps/30to5min/';
if(exist(outputDir,'dir') == 0)
	mkdir(outputDir);
end

% Open video object
vidObj = VideoWriter(sprintf('%s%s-%s-%s-%s', outputDir, date(5:6), date(7:8), date(1:4), startstr), 'MPEG-4');
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);

forecast = load(sprintf('%s/%s',inputDir,files(file_id).name));

for i = file_id:file_id+9
	nowcast = load(sprintf('%s/%s',inputDir,files(i+1).name));
	figure
	% Plot nowcast as solid base
	imagesc(nowcast.advect{1}, [-3 6]);
	colormap(ovcmap)
	hold on
	
	% Plot forecast as overlay -- set AlphaData to make everything
	% transparent except clouds
	
	% Search for index corresponding to nowcast in forecast struct
	q = find(forecast.time == nowcast.time(1));
	
	forecastadvect = imagesc(forecast.advect{q}.*3, [-3 6]);
	ovalpha = get(forecastadvect, 'CData');
	ovalpha(ovalpha <= 0 | ovalpha == 4 | ovalpha == 5) = 0;
	ovalpha(ovalpha == 1 | ovalpha == 2 | ovalpha == 3 | ovalpha == 6) = 0.4;
	colormap(ovcmap)
	set(forecastadvect, 'AlphaData', ovalpha)
	axis image
	axis off
	title(sprintf('%is Forecast Cloud Map overlaid on Nowcast Cloud Map [%s (UTC)]', (q-1)*30, datestr(nowcast.time(1))))
	hold off
	set(gcf,'Position',[520 210 640 620])
	currFrame=getframe(gcf);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	writeVideo(vidObj,currFrame);
	close
end
close(vidObj)
end	