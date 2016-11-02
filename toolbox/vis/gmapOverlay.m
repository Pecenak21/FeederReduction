function gmapOverlay(img, lonlim, latlim, alpha_, option)
% overlay an image on top of the google map with specified transparency
% 
% input:
%			img		
%			latlim	: latitude limit matrix in degree. E.g. latlim = [37.8 37.9]
%			lonlim	: latitude limit matrix in degree. E.g. lonlim = [-117.2534 -117.2127]
%			alpha_	: (optional) transparency. Default: .7
%
% option setting:
%		option.type : support option.type = 'cloud'
%
% Example:
%			op.type = 'cloud';
%			gmapOverlay(img, [-117.2534 -117.2127], [32.8695 32.8928], .7, op);

if ~exist('alpha_','var') || isempty(alpha_)
	alpha_ = 0.7;
end

if ~exist('option','var') || isempty(option) || ~isfield(option, 'type')
	option.type = 'cloud';
end

% georeference
ucsd = siDeployment('UCSD');
% georeference
R = genGeoRef(ucsd.ground);
% R = georasterref;
[a, b, c] = size(img);
% R.RasterSize = size([a b]);
% R.Latlim = latlim;
% R.Lonlim = lonlim;

figure;
set(gcf,'visible','on')

% plot google map
plot(R.Lonlim,R.Latlim,'.r','MarkerSize',1)
plot_google_map('MapType','hybrid')
hold on

set(gcf,'position',[0 0 1200 900]);
set(gcf,'color','w')
grid on; box on;
set(gca,'fontsize',20); xlabel('Latitude');ylabel('Longitude');

% initilize data to plot
x = zeros(a,b);
y = zeros(a, b, 3 );

switch option.type
	case 'cloud'
		% cloud
		x(img>0) = 1;
		% undefined region
		x(img<0) = nan;
		% blue sky
		x2 = x;
		x2(x2==0) = alpha_;
		y(:,:,1) = x;y(:,:,2) = x2;y(:,:,3) = x2;
	otherwise
		x = img(:,:,1);
		y = img;
end
	
almap = zeros(a,b);
almap(x>0) = .7;
mapshow(y,R,'AlphaData',almap);
xlim([lonlim(1)-.0001,lonlim(2)+.0001]);
ylim([latlim(1)-.0001,latlim(2)+.0001]);

end