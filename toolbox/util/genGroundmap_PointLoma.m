%% Generate groundmap for Redlands
% use USI_1_6 as reference location
u6 = siImager('USI_1_9');

%% get initial default deployment file
de = siDeployment('PointLoma_UCSD');

%% load footprint excel file which has the whole site footprint and PV footprint
load('data/ground.mat');

%% distance from usi_1_6 to generate the groundmap
% Go North (round up to 100m unit)
d(1) = domain.n;
% Go East
d(2) = domain.e;
% Go South
d(3) = domain.s;
% Go West
d(4) = domain.w;
%% resolution
res = 2.5; % in m

%% generate footprint
dx = size(ground.latitude,2);
dy = size(ground.latitude,1);

% irradiance sensor (GHI)
limx = [min(ground.longitude(:)) max(ground.longitude(:))];
limy = [min(ground.latitude(:)) max(ground.latitude(:))];

%% PV sites
footprint = load('data/footprint.mat');
count = 0;
pv=pvsystem;
for i = 1:length(footprint.pv)
	count = count + 1;
	pv(count) = pvsystem;
	pv(count).cenlat = footprint.pv(i).pos.latitude;
	pv(count).cenlon = footprint.pv(i).pos.longitude;
    pv(count).kVA = footprint.pv(i).scaleFactor;
	pv(count).calArea;
	% name
	% PV location
	pixelx = round( dx * (pv(count).cenlon - limx(1)) / (limx(2)-limx(1)) );
	pixely = round( dy * (pv(count).cenlat - limy(1)) / (limy(2)-limy(1)) );
	
	idx =	(ground.latitude > pv(count).lat(1)) & (ground.latitude < pv(count).lat(2)) & ...
			(ground.longitude > pv(count).lon(1)) & (ground.longitude < pv(count).lon(2)) ;
	if sum(idx(:)) > 0
		footprint.inverter(idx) = i;
	else
		footprint.inverter(pixely,pixelx) = i;
	end
end
footprint.pv = pv;

save('data/footprinttest.mat','-struct','footprint');