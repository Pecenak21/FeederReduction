%% This generates ground map
% IMPORTANT! I RAN THE FIRST BLOCK OF YOUR PointLoma_test.m IN ORDER TO GET lat AND lon ARRAYS. IT CONTAINED YOUR SDGE FEEDER DATA.

%%%%% Run originally to generate ground map. Change parameters as needed
%clear all
position.latitude = 32.7251;
position.longitude = -117.2488;
position.altitude = 93; % [m]
d = [2500 1000 2500 1000]; % domain size, centered on position struct [N E S W]
pl = geo_getGroundmap( position, d, 2.5); % generate ground map
save('PointLoma_forecast/PL_ground.mat', '-struct', 'pl')
%clear all
%%%%%

%%%%% Load previously generated ground map instead of generating again

target.ground = load('PointLoma_forecast/PL_ground.mat'); % Stick PL_ground.mat into deployment folder, and rename it to 'ground.mat'

% Name this the same as your data type for power reference (GHI or PV power output). I used "none" because I did not have any PV data.
% Whatever you name this, name it the same thing as what you put in the data_type field in the deployment.conf
% E.g.: deployment.conf:        data_type              AC
% AC = NaN(...);
none = NaN(size(target.ground.longitude)); % Blank footprint variable

%% load feeder data to get lat lon information of the loads
f = load('PointLoma_forecast/SDGE_feeders.mat','f480');f = f.f480;

% find max lon, lat
lon = [f.X]; lat = [f.Y];

%% TODO Maxime: make footprint of PV panels
% load PV locations

% put it into the grid with a specific id for each PV system 

%% Make footprint
ulat = unique(lat(:),'stable');
ulon = unique(lon(:),'stable');

for id = 1:length(ulat)
	% Label in footprint
	[r(id), c(id)] = geo_getLatLonIndex(ulon(id), ulat(id), target.ground.longitude, target.ground.latitude); % Map lat/lon to row and column within footprint
	% IMPORTANT: CHECK YOUR OUTPUT AFTER. IF YOU GET AN r OR c AT THE MINIMUM/MAXIMUM SIZE OF YOUR DOMAIN, YOUR DOMAIN IS TOO SMALL.
	
	none(r(id), c(id)) = 100 + id;	% IMPORTANT: Convention within forecast code is for PV panels to begin from 101. Also, should change variable name here.
									% Should just change "none" above, and shift+enter it.
	
	% PV position
	% Everything is default value for now, but you are able to change each pv panel individually. I'll leave how to do that up to you.
	pv(id).pos.longitude = ulon(id);
	pv(id).pos.latitude = ulat(id);
	pv(id).pos.altitude = 100; % Default value
	pv(id).tilt = 0; % Default
	pv(id).azimuth = 180; % Default
	
	% Power model. Currently models EVERY PV system the same way. I'll let you figure out how to set them differently individually.
	pv(id).eff_coeff = [-4.6986e-18, 3.7089e-13, -1.1256e-08, 1.4799e-04, 0.3500]; % Polynomial power efficiency curve fit, generated using pvfit.m
	pv(id).scaleFactor = 22.1960; % Approximately in kW
	
	% TODO: replace this by PV name
	inverterNames{id} = num2str(id); 
end

save('PointLoma_forecast/PL_footprint.mat', 'none', 'pv', 'inverterNames') % Tries to save your stuff into a .mat file, obviously. Be sure to change 'none' to whatever your footprint field was.
% Stick PL_footprint.mat into the deployment folder, and rename it 'footprint.mat'
% I think I also saved the r and c variables as 'XY.mat' for plotting purposes. You may find this helpful if you want to use prettyplot.