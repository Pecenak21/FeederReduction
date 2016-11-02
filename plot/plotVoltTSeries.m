function vpath = plotVoltTSeries(res,prefix)
% plot voltage profile for the whole day and make movie out of the resulted images
% input: 
%       res:        result structure from simMain or path to file that contain the result struct
%       prefix:     (optional) prefix to saved filenames

if ischar(res), res=load(res); end
%if ~exist('outDir','var'), outDir = 'plotVolt'; end
%if ~exist(outDir,'dir'), mkdir(outDir); end
if ~exist('prefix','var'), prefix = 'voltProf'; end

t = res.time;
% generate all 
yl = [floor(min(res.Voltage(:))*100)/100 ceil(max(res.Voltage(:))*100)/100];
pf = fNamePrefix(prefix,res); pf = pf{1};
if ~exist(pf,'dir'), mkdir(pf); end
for i = 1:length(t)
    figure(100); 
    plot(res.Dist,res.Voltage(i,:),'.');
    xlabel('Distance from substation [km]');
    ylabel('Voltage [p.u.]');
    title(datestr(t(i),'yyyy-mm-dd HH:MM:SS'));
    ylim(yl); grid on; box on;
    saveas(gcf,[pf '/' sprintf('%04d',i) '.png']);
end

% make video
ffmpegtranscode([pf '/%04d.png'],[pf '.mp4'],'InputFrameRate',60,...
         'x264Tune','animation');
vpath = [pf '.mp4'];
end