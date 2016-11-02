%% Generate shadow map -google map overlay video on Point Loma feeder - 480 Cabrillio
% forecast data dir
% def_addpath
% datadir = '/mnt/lab_18tb3/database/USI/analysis/LateJanRevision';
datadir = 'P:/database/USI/analysis/LateJanRevision';

% load feeder data
f = load('data/SDGE_feeders.mat','f480');f = f.f480;

% find max lon, lat
lon = [f.X]; lat = [f.Y];
% grid limit
g.ylim = [min(lat(:)) - .005, max(lat(:)) + .005 ];
g.xlim = [min(lon(:)) - .02, max(lon(:)) + .02 ];

% UCSD imager
u2 = siImager('USI_1_2');
% UCSD ground
ucsd = siDeployment('UCSD');
% georeference
R = genGeoRef(ucsd.ground);
R.Latlim = g.ylim;
R.Lonlim = g.xlim;

filename = 'tmp_fig/shadowOnPointLomaFeeder480_2.gif';

%%
% for t = datenum([2012 12 14 18 00 00]):30/3600/24:datenum([2012 12 14 20 00 00])
for t = datenum([2012 11 14 20 00 00]):30/3600/24:datenum([2012 11 14 21 00 00])
	tic
	%%
	figure(1)
	set(gcf,'visible','on','color','w')
	xlim(g.xlim);
 	ylim(g.ylim); 
	
	% plot google map
	plot(R.Lonlim,R.Latlim,'.r','MarkerSize',1)
	plot_google_map('MapType','hybrid')
	hold on
	
	% set plot format
	xlim(g.xlim);
 	ylim(g.ylim); 
	set(gcf,'position',[0 0 900 900]);
	grid on; box on;
	set(gca,'fontsize',20); xlabel('Latitude');ylabel('Longitude');
	set(gca,'xtick',[g.xlim(1) : (g.xlim(2)-g.xlim(1)) /3 : g.xlim(2)]);
	set(gca,'ytick',[g.ylim(1) : (g.ylim(2)-g.ylim(1)) /5 : g.ylim(2)]);
	
	% plot feeder map
	for i = 1:length(f)
		plot([f(i).X]',[f(i).Y]','g','linewidth',2);
	end
	
	% load data
% 	fc = loadforecastdata(u2.t,t,'forecast',
	fc = load(['P:\database\USI\analysis\LateMay2013_v2\20121114/forecast/forecast_' datestr(t,'yyyymmddHHMMSS') '.mat'],'shadow');
	d = fc.shadow{1};
	d = medfilt2(d, [7 7]);
	ch = load(['P:\database\USI\analysis\LateMay2013_v2\20121114/cloudheight/cloudheight_' datestr(t,'yyyymmddHHMMSS') '.mat']);
	height = ch.height;
	% overlay data
		
	% initilize data to plot
	[a b] = size(d);
	clear x y;
	x = zeros(a,b);
	y = zeros(  a, b, 3 );
	
	% cloud
	x(d>0) = 1;
	% undefined region
	x(d<0) = nan;
	% blue sky
	x2 = x;
	x2(x2==0) = .7;

	y(:,:,1) = x;y(:,:,2) = x2;y(:,:,3) = x2;
% 	almap = uint8(255*ones(a,b));
	almap = zeros(a,b);
 	almap(d==0) =0;
	almap(d>0) = .7;
 	mapshow(y,R,'AlphaData',almap); 
	xlim([-117.2766 -117.2209]);
 	ylim([32.7024   32.7478]); 
	grid on; box on;
	set(gca,'fontsize',20); xlabel('Latitude');ylabel('Longitude');
	set(gca,'xtick',[round(100*g.xlim(1))/100 : round(100*(g.xlim(2)-g.xlim(1)) /5)/100 : round(g.xlim(2)*100)/100]);
	set(gca,'ytick',[round(100*g.ylim(1))/100 : round(100*(g.ylim(2)-g.ylim(1)) /7)/100 : round(100*g.ylim(2))/100]);
	
	%% save image
	text(0.05,.1, sprintf('Cloud Height: %.0f m',height),'units','normalized','Fontsize',20,'color','r');
	text(0.05,.05, sprintf('Time (PST): %s ',datestr(t-8/24,'HH:MM:SS mmm-dd-yyyy')),'units','normalized','Fontsize',20,'color','r');
% 	saveas(gcf,sprintf('tmp_fig/decmapbinary_normalized_%.0fm.png',height))
	
	% make video
	drawnow
	frame = getframe(1);
	im = frame2im(frame);
	[imind,cm] = rgb2ind(im,256);
	if t == datenum([2012 12 14 18 00 00])
	imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
	else
	imwrite(imind,cm,filename,'gif','WriteMode','append');
	end
	hold off
	toc
% 	pause
end