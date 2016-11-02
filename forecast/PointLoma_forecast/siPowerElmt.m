%% Solar Resource Assessment
%  Power plant analysis
%
%  Title: Power Forecasting Element
%
%  Description:
%    This function will extract a set time series of data from the power plant
%    as indicated by the time input, and it will then perform conversions of the
%    data into horizontal clear sky index using the Muneer and Boland models.
%    This data will then be returned in a uniform structure.
%
%  Input:
%	 target		siDeployment struct
%	 time		times for which to get power
%				Specify a long list to get power every 1 sec between the last two times in the list
%				Specify a single time to get power every 1 sec for the 30 seconds leading up to it
%	 scaleFactor (optional) to convert readings in kW to W or vice-versa if desired
%				TODO: scaleFactor should probably be a property of the deployment setup?
%
function pwr_elmt = siPowerElmt( target, time, scaleFactor )
%% Process Input Arguments

% Note: DEMROES data is used and it's in -8 hour timezone (LA)
conf.timezoneloc = 'Los Angeles';
conf.timezoneDST = 0;

% lookback time
lookback = 30/3600; % in hours

% Data resolution (in seconds). This is possibly a deployment setting.
dres = 1; %[s]

% Image time resolution default
ires = 30; %[s]

% Default is to convert kW readings to W.
if( nargin < 3 || isempty(scaleFactor) )
	scaleFactor = 1000; %[W/kW]
end
% if no time delta is specified, use 30 seconds
if( length(time) < 2)
	time = time - [ ires/(24*3600) , 0 ];
end

%% Store the time
pwr_elmt.time = (time(end) - ires/24/3600 + dres/(24*3600)):dres/(24*3600):time(end);

%% Retrieve Data
if (strcmpi(target.name,'ucsd') || strcmpi(target.name,'ucsd_augmented')) && isfield(target,'data_type') && strcmpi(target.data_type(1),'g')
	if strcmpi(target.data_type, 'ghi'), fprintf(' Ground data type: GHI sensor\n'), end
	if strcmpi(target.data_type, 'ghipv'), fprintf(' Ground data type: GHI sensor & PV\n'), end

	%% Do DEMROES sensors
	% Fetch DEMROES data, compute clear sky model GHI, divide measured GHI by modeled clear sky GHI to obtain clear sky index kt, to be used as
	% input into siForecastGHI's kt PDF procedure.
	[pwr_elmt.sensor.GHI, pwr_elmt.time] = getDEMROES( (time(end) - lookback/24 + 1/24/3600), time(end), target);
	csk = clearSkyIrradiance( target.ground.position , pwr_elmt.time, target.tilt, target.azimuth );
	pwr_elmt.sensor.kt = pwr_elmt.sensor.GHI  ./ repmat( csk.gi(:) , [1 size(pwr_elmt.sensor.GHI,2)] );
	
	%% Do PV arrays if present
	% Currently assumes PV arrays located on UCSD on one PI server. Fetch PV data from UCSD PI server, compute plane of array clear sky
	% model GI. Model clear sky PV power using modeled clear sky GI and empirical efficiency curve (currently 4th order polynomial fit).
	% Divide measured PV power by modeled PV power to obtain kt.
	if strcmpi(target.data_type,'ghipv')
		pwr_elmt.inverter.power = target.ucsdpi.getPower(time(end-1)+dres/(24*3600), time(end),1);
		csk.pv.gi = ghi2gi( target.footprint.pv.azimuth , target.footprint.pv.tilt , target.footprint.pv.pos , csk.time , csk.ghi );
		clrPwrModel = csk.pv.gi'*target.footprint.pv.scaleFactor;
		clrPwrModel = polyval(target.footprint.pv.eff_coeff, clrPwrModel).*clrPwrModel;
		pwr_elmt.inverter.kt = pwr_elmt.inverter.power./clrPwrModel;
	end
elseif strcmpi(target.name,'redlands')
	%% we have 1 GHI data from a GHI sensor at PV022 and 4 power data sets. We will convert 4 power data sets to GHI and use them as GHI sensors
	%% load GHI data from the irradiance sensor
	
	%% convert GHI to kt
	
	%% load power data from the 4 power data sets
	% check if data exist
	% convert to kt
	
	%% concat kt for all 5 stations
	pwr_elmt.time2h = ( pwr_elmt.time(end) - lookback/24 + 1/24/3600 ) :1/24/3600: pwr_elmt.time(end);
	pwr_elmt.sensor.kt;
	
elseif isfield(target,'data_type') && strcmpi(target.data_type(1),'n') % No data
	pwr_elmt.inverter.kt = nan(30,1); % Feed in 30 seconds worth of NaN. Default kt values will be used.
	
else  %% works for Henderson
	% Get all the data between the last time step and this time step. Here the
	% presumed data resolution is 1sec.
	inverter = target.pi.obj.getPower(time(end-1)+dres/(24*3600), time(end),1);

	% Convert to matrix, time by inverter. Scale to SI units.
	pwr_elmt.inverter.power = cat(2,inverter(:).power)' * scaleFactor;

	%% Compute gi
	pwr_elmt.inverter.gi  = pwrPowerToGI( pwr_elmt.inverter.power' , target.design.power )';

	%% Compute the plane of array kt
	% need the clear sky irradiance first
	csk = clearSkyIrradiance( target.ground.position , pwr_elmt.time, target.tilt, target.azimuth );
	% dividing by clear sky gives kt
	pwr_elmt.inverter.kt = pwr_elmt.inverter.gi  ./ repmat( csk.gi , [1 size(pwr_elmt.inverter.gi,2)] );
end

end
