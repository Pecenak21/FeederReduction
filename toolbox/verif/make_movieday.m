% Makes a 30 fps, 525 x 525 .avi of high exposure raw images. Syntax: make_movie(YYYYMMDD, USI#, 'HH (UTC)')
% E-mail Handa if you have any questions about this code

function make_movieday(USI_number, day)

outputDir = 'Z:/infobase/daily_video/';
if(exist(outputDir,'dir') == 0)
	mkdir(outputDir);
end

if ~ischar(day)
	day = num2str(day);
end

vidObj = VideoWriter(sprintf('%susi1-%i/%s-%s-%s.avi', outputDir, USI_number, day(1:4), day(5:6), day(7:8)));
vidObj.FrameRate=30;
open(vidObj);

for hour_start = 15:23
	
	fprintf('Date: %s/%s/%s, USI Number: %i, Hour (UTC): %i\n', day(5:6), day(7:8), day(1:4), USI_number, hour_start)
	b=dir(sprintf('Z:/database/USI/images/usi1-%i/%s/%s%i*_prev.jpg', USI_number, day, day, hour_start));

	for i=1:3:length(b)

		I=imread(sprintf('Z:/database/USI/images/usi1-%i/%s/%s', USI_number, day, b(i).name));

		imshow(I,'InitialMagnification',60)
		text(10,size(I,2)/2,'E','FontSize',18,'Color','red')
		text(size(I,1)-30,size(I,2)/2,'W','FontSize',18,'Color','red')
		text(size(I,1)/2,30,'N','FontSize',18,'Color','red')
		text(size(I,1)/2,size(I,2)-30,'S','FontSize',18,'Color','red')
		text(size(I,2)*4/6+160,size(I,2)-20,[sprintf('%02d',str2num(b(i).name(9:10))),':',b(i).name(11:12),':',b(i).name(13:14)],'FontSize',14,'Color','green')
		text(10,size(I,2)-20,[day(5:6),'/',day(7:8),'/',day(1:4)],'FontSize',14,'Color','green')
		currFrame=getframe;
		writeVideo(vidObj,currFrame);
	end
end

%Close file
close(vidObj)
fprintf('Creation of Z:/infobase/daily_video/usi1-%i/%s-%s-%s.avi completed.\n', USI_number, day(1:4), day(5:6), day(7:8))


% ______________________________________
% 
% But you will have to change the folder names, etc. The most important
% lines if you want to make your own version are:
% 
% %Open the video object with a given name and frame rate
% 
% vidObj = VideoWriter(['name.avi']);
% vidObj.FrameRate=30;
% open(vidObj);
% 
% % show the USI image
%     imshow(I)
% 
% %Capture the frame and write it on the movie
%     currFrame=getframe;
%     writeVideo(vidObj,currFrame);
% 
% %Close the video object
% 
% close(vidObj)