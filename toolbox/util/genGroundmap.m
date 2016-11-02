%% Generate groundmap for Redlands
% use USI_1_6 as reference location
u6 = siImager('USI_1_6');

%% get initial default deployment file
de = siDeployment('Redlands');

%% load footprint excel file which has the whole site footprint and PV footprint
fp = excel2obj([de.source '/RedlandsSiteFootprint.xlsx']); fp = fp.FootPrint;

%% distance from usi_1_6 to generate the groundmap
% get whole site lat, lon
[id, id]= ismember('WholeSite',{fp.Site}');
% Go North (round up to 100m unit)
d(1) = 1000*greatCircleDistance(u6.position.latitude,u6.position.longitude,fp(id).MaxLat,u6.position.longitude);
d(1) = ceil(d(1)/100)*100;
% Go East
d(2) = 1000*greatCircleDistance(u6.position.latitude,u6.position.longitude,u6.position.latitude,fp(id).MaxLon);
d(2) = ceil(d(2)/100)*100;
% Go South
d(3) = 1000*greatCircleDistance(u6.position.latitude,u6.position.longitude,fp(id).MinLat,u6.position.longitude);
d(3) = ceil(d(3)/100)*100;
% Go West
d(4) = 1000*greatCircleDistance(u6.position.latitude,u6.position.longitude,u6.position.latitude,fp(id).MinLon);
d(4) = ceil(d(4)/100)*100;

%% resolution
res = 2.5; % in m

%% generate groundmap 
gm = geo_getGroundmap(u6.position,d,res);
gm.irradianceSensor.lat = fp(i).Lat_Center_;
gm.irradianceSensor.lon = fp(i).Lon_Center_;
gm.irradianceSensor.name = fp(i).Site;
gm.irradianceSensor.address = fp(i).Address;

save([de.source '/groundmap.mat'],'-struct','gm');

%% generate footprint
dx = size(gm.latitude,2);
dy = size(gm.latitude,1);

% irradiance sensor (GHI)
[id, id] = ismember('IrradianceSensor',{fp.Site}');
footprint.GHI = zeros(dy,dx);
footprint.GHInames = {'SPVP022_IrrSensor_GHI'};
limx = [min(gm.longitude(:)) max(gm.longitude(:))];
limy = [min(gm.latitude(:)) max(gm.latitude(:))];
pixelx = round( dx * (fp(id).Lon_Center_ - limx(1)) / (limx(2)-limx(1)) );
pixely = round( dy * (fp(id).Lat_Center_ - limy(1)) / (limy(2)-limy(1)) );
footprint.GHI(pixely,pixelx) = 1;

%% PV sites
footprint.inverter = zeros(dy,dx);
id = find(~cellfun(@isempty,regexpi({fp.Site}','SPVP.*')));
count = 0;
for i = id'
	count = count + 1;
	pv(count) = pvsystem();
	pv(count).name = fp(i).Site;
	pv(count).cenlat = fp(i).Lat_Center_;
	pv(count).cenlon = fp(i).Lon_Center_;
	pv(count).lat = [fp(i).MinLat fp(i).MaxLat];
	pv(count).lon = [fp(i).MinLon fp(i).MaxLon];
	pv(count).kVA = fp(i).AC_Rating_MW_;
	pv(count).calArea;
	pv(count).note = fp(i).Address;
	% name
	footprint.inverternames{count} = fp(i).Site;
	% PV location
	pixelx = round( dx * (pv(count).cenlon - limx(1)) / (limx(2)-limx(1)) );
	pixely = round( dy * (pv(count).cenlat - limy(1)) / (limy(2)-limy(1)) );
	
	idx =	(gm.latitude > pv(count).lat(1)) & (gm.latitude < pv(count).lat(2)) & ...
			(gm.longitude > pv(count).lon(1)) & (gm.longitude < pv(count).lon(2)) ;
	if sum(idx(:)) > 0
		footprint.inverter(idx) = i;
	else
		footprint.inverter(pixely,pixelx) = i;
	end
end
footprint.pv = pv;

save([de.source '/footprint.mat'],'-struct','footprint');