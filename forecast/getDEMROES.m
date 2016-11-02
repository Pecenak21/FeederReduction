function [GHI, timeUTC, stationInfor, GHIavg, GHI_kt_avg, persistenceGHI, GHIclrsky] = getDEMROES(starttime, endtime, target, station, horizon)
% Get DEMROES GHI data
%	DEMROES is the UCSD weather sensor network.
%
% Inputs:
%			starttime : either time string or datenum in UTC
%			endtime : either time string or datenum in UTC
%			target : siDeployment structure for UCSD (optional, defaults to 'UCSD')
%			station (optional): either station string name or numerical id. Default: (all stations).
%			
% Output:
%			GHI
%			timeUTC
%			stationInfor : names of the stations
%			GHIavg
%			GHI_kt	: of GHIavg
%
% Example:
%			stime = datenum([2012 12 15 08 00 00]);
%			etime = datenum('2012-12-16 08:30:00');
%			ghi = getDEMROES(stime, etime); % get all stations
%			ghi = getDEMROES(stime, etime, [], 'MOCC'); % for MOCC
%			ghi = getDEMROES(stime, etime, [], [1 2 3]); % for fisrt 3 stations
%
% UCSD DEMROES stations: 'BMSB','CMRR','EBU2','HUBB','MOCC','POSL'
% UCSD augmented stations:  'TIOG', 'MWFM', 'VAMC', 'FRST', 'CANY'

if ischar(starttime)
	starttime = datenum(starttime);
end

if ischar(endtime)
	endtime = datenum(endtime);
end

if ~exist('target','var') || isempty(target)
	% Assume a default target of UCSD.  This is for quick hacking purposes to making the function easy to call.
	% If you are writing a script, you should generally actually get an appropriate deployment and use it rather than relying on the default.  Especially because that will re-load it from disk every time.
	warning('getDEMROES:missingDeployment','Please input a deployment, e.g. target = siDeployment(''UCSD'')');
	target = siDeployment('UCSD');
end
if ~exist('horizon','var'), horizon = 0; end

% Define site info

site = struct('name',target.footprint.GHInames,'timezoneloc',target.footprint.timezoneloc,'timezoneDST',num2cell(target.footprint.timezoneDST),'GHIcal',num2cell(target.footprint.GHIcal));
% This code could allow us to do something like this automatically, but
% it's not needed in general (only need MOCC), and gives different results
% than what Handa was originally using to calculate shading at MOCC, so for
% now, we'll just do MOCC:
% for sid = 1:length(site)
%     m = (target.footprint.GHI==sid);
%     site(sid).longitude = mean(target.ground.longitude(m));
%     site(sid).latitude = mean(target.ground.latitude(m));
%     site(sid).altitude = mean(target.ground.altitude(m));
% end
% Hard-code MOCC position: currently only one used. If necessary, positions can be added to deployment.
id = find(strcmp('MOCC',{site.name}));
if(~isempty(id))
    site(id).longitude = -117.222547;
    site(id).latitude = 32.878443;
    site(id).altitude = 103;
end

if ~exist('station','var') || isempty(station)
	sid = 1:length(site);
else
	if ischar(station) || iscellstr(station)
		[~,sid] = ismember(lower(station),lower({site.name}));
		if(any(~sid))
			error('getDEMROES:unknown_station','Unknown station: %s\n',station{~sid});
		end
	else
		sid = station; % num
	end
end
stationInfor = site(sid);

% Define 1s DEMROES data path
setup = readConf(getConfPath('setup'),0);
databasePath = normalizePath(setup.DEMROES_DATA_DIR);

% DST fix, see below
% start for most stations, end for most stations, start for CMRR, start for TIOG
dstfixt = datenum([2014,04,18, 01,00,00;  2014,07,28, 08,00,00;  2014,04,28, 21,30,00;  2014,04,18, 03,00,00]);
dstfixtimes = repmat(dstfixt(1:2)',numel(site),1);
% start, end, cmrrstart, tiogstart
posl_id = find(strcmp('POSL',{site.name})); if(isempty(posl_id)), posl_id = nan; end
dstfixtimes(strcmp('CMRR',{site.name}),1) = dstfixt(3);
dstfixtimes(strcmp('TIOG',{site.name}),1) = dstfixt(4);
% except  for CMRR

%% Load DEMROES GHI data
% Initialize
t = starttime - max(horizon)/24/60 : 1/3600/24 : endtime; % time stamp with increment of 1 second; Read in GHI data 15 minutes prior (for error plot)
GHI = nan( length(t), length(site) );
for id = sid
	if(endtime > nowUTC - 5/60/24) %If the endtime is within 5 minutes of now, we can safely assume the operational forecast is in effect
		% Setup for Demroes database
		dbName = 'demroes';
		dbUser = 'matlabuser';
		dbPassword = 'matlabuser_psswrd';
		dbDriver = 'com.mysql.jdbc.Driver';
		dbUrl = 'jdbc:mysql://panther8.ucsd.edu/demroes';
		% Open the connection
		try
			usedbtoolbox = true;
			conn = database(dbName,dbUser,dbPassword, dbDriver,dbUrl);
		catch err
			if(~strcmp(err.identifier, 'MATLAB:license:checkouterror')), rethrow(err); end
			usedbtoolbox = false;
			conn = com.mathworks.toolbox.database.databaseConnect(dbName,dbUser,dbPassword, dbDriver,dbUrl);
			v = conn.makeDatabaseConnection();
			if(v.get(0)) % worked
				conn = conn.getValidConnection(v);
			else
				error('unable to connect to DB server');
			end
		end

		% Setup for the site and the period under question
		%queryString = sprintf('select TmStamp, %s from %s_HighFreq where TmStamp > ''%s'' AND TmStamp <= ''%s''', colName, site(id).name, datestr(toLocalTime(starttime,site(id).timezoneloc,site(id).timezoneDST),31), datestr(toLocalTime(endtime,site(id).timezoneloc,site(id).timezoneDST),31));
		% Can get Unix timestamp like this, and then convert as datenum('1970-01-01 00:00:00') + GHIdata.Unixtime/24/3600;
		queryString = sprintf('select UNIX_TIMESTAMP(TmStamp), Solar_Radiation from %s_HighFreq where TmStamp > ''%s'' AND TmStamp <= ''%s''', site(id).name, datestr(toLocalTime(starttime,site(id).timezoneloc,site(id).timezoneDST),31), datestr(toLocalTime(endtime,site(id).timezoneloc,site(id).timezoneDST),31));
		% Query the data
		if(usedbtoolbox)
			curs = exec(conn, queryString);
			curs = fetch(curs);
			dbdata = curs.Data; % is a cell array with the data you wanted in it.
		else
			% curs = com.mathworks.toolbox.database.sqlExec(queryString, conn);
			% fet = com.mathworks.toolbox.database.fetchTheData(conn, curs, queryString);
			curs = conn.prepareCall(queryString);
			dbdata = curs.executeQuery();
		end
		
		% Re-formatting the data to match USI conventions
		if(~usedbtoolbox && ~dbdata.wasNull)
			dbdata.last(); n = dbdata.getRow();
			dbdataraw = nan(n,2);
			dbdata.beforeFirst();
			i = 0;
			while(dbdata.next())
				i = i+1;
				dbdataraw(i,1) = dbdata.getLong(1);
				dbdataraw(i,2) = dbdata.getDouble(2);
			end
			if(i~=n)
				warning('number of rows fetched does not match number of rows advertised');
			end
			GHIdata.time_day = toUTC( datenum([1970 1 1 0 0 0]) + dbdataraw(:,1)/24/3600, site(id).timezoneloc, site(id).timezoneDST );
			GHIdata.GHI_day = dbdataraw(:,2);
		elseif(usedbtoolbox && iscell(dbdata) && ~strcmp(dbdata{1},'No Data'))
			%GHIdata.time_day = toUTC( datenum(vertcat(dbdata{:,1}),'yyyy-mm-dd HH:MM:SS'), site(id).timezoneloc, site(id).timezoneDST );
			GHIdata.time_day = toUTC( datenum([1970 1 1 0 0 0]) + vertcat(dbdata{:,1})/24/3600, site(id).timezoneloc, site(id).timezoneDST );
			GHIdata.GHI_day = vertcat(dbdata{:,2});
		else
			GHIdata.time_day = [];
			GHIdata.GHI_day = []; % No data found
		end
	else % Otherwise continue with acquiring Demroes for the historical forecast
		% change time to DST format (DEMROES format) to find the days and get correct data (data files are named based on DST format aaaaaaa)
		day = floor(toLocalTime(starttime,site(id).timezoneloc,site(id).timezoneDST)):...
			floor(toLocalTime(endtime,site(id).timezoneloc,site(id).timezoneDST));
		% calculate day of year by days since the previous year's Dec 31
		doy = datevec(day(:));
		doy = day - datenum([doy(:,1)-1, repmat([12, 31],size(doy,1),1)])';
		GHIdata = struct('time_day',[],'GHI_day',[]);
		for dn = 1:length(day)
			d = day(dn); doy_ = doy(dn);
			fn = sprintf('%s/%s/%s_%s_%i.mat', databasePath, site(id).name, site(id).name, datestr(d, 'yyyy'), doy_);
			if exist(fn,'file')
				% Load GHI data
				GHIdata_ = load(fn);
				GHIdata.time_day = vertcat(GHIdata.time_day, toUTC( GHIdata_.time_day(:) , site(id).timezoneloc, site(id).timezoneDST ));
				if isfield(GHIdata_, 'GHI_day') % DEMROES convention
					GHIdata.GHI_day = vertcat(GHIdata.GHI_day, GHIdata_.GHI_day(:));
				elseif isfield(GHIdata, 'data_day') % Stealth pyranometer convention
					GHIdata.GHI_day = vertcat(GHIdata.GHI_day, GHIdata_.data_day(:));
				else
					error('getDEMROES:missingField', 'Error in file %s. Cannot find field "GHI_day" or "data_day".', fn)
				end
			end
		end
	end
	% Filter MOCC data by SZA only before Jan 18
	if( (id == 5) && (starttime < datenum([2013 01 18 08 00 00])) )
		MOCC_SUN = siSunPosition( GHIdata.time_day , site(id) ); % all that matters is that we use a struct with latitude/longitude/altitude fields
		MOCC_maxzenith = atand(1.41 ./ cosd( abs(180 - [MOCC_SUN.azimuth]) ) );

		% only filter data before Jan 18
		GHIdata.GHI_day( (MOCC_SUN.zenith(:) > MOCC_maxzenith) & (GHIdata.time_day < datenum([2013 01 18 08 00 00])) ) = NaN;
	end
	% End MOCC filter
	
	if(isempty(GHIdata.time_day)), continue; end
	
	% Timestamps were incorrectly recorded in PDT instead of PST from 2014-04-17 17:35 to 2014-07-28 08:00.  Fix that:
	if (id~=posl_id && GHIdata.time_day(end) > dstfixtimes(id,1) && GHIdata.time_day(1) < dstfixtimes(id,2))
		mask = GHIdata.time_day > dstfixtimes(id,1) & GHIdata.time_day < dstfixtimes(id,2);
		GHIdata.time_day(mask) = GHIdata.time_day(mask) - 1/24;
	end

	% find matched time
	[lia,lib] = ismember( round(GHIdata.time_day.*(24*3600)), round( t.*(24*3600)) );
	if sum(lia) > 0
		lib(lib==0) = [];
		% assign corresponding GHI
		GHI(lib,id) = GHIdata.GHI_day( lia );
	end
	GHI(:,id) = GHI(:,id)*site(id).GHIcal; % Multiply pyranometer data by calibration factor
end

GHI = GHI(:,sid);

% clean up GHI data by removing all timestamps that does not have any data from any station
id = ~isnan(GHI);
id = sum(id,2) > 0;

timeUTC = t(id);
GHI = GHI(id,:);
% Average GHI data if needed
if sum(id) > 0
	% refine one more time based on clear sky irradiance
	csk = clearSkyIrradiance( target.ground.position , timeUTC, target.tilt, target.azimuth );
	id = csk.gi > 0;
	timeUTC = timeUTC(id);
	GHI = GHI(id,:);
	% Take average, change nan values to 0, if asked for (i.e. averages are requested)
	if nargout > 3
		GHIavg = nanmean(GHI,2);
		% Convert DEMROES GHI Data to kt
		GHI_kt_avg = GHIavg./csk.gi(id)';
		GHIclrsky = csk.gi(id);
		
		% apply a bit of quality control here
		GHI_kt_avg(GHI_kt_avg>1.5) = 1.5;
	end
	% Compute persistence values if requested (mostly used for prettyplot)
	if nargout > 5
		persistenceGHI.GHI_kt_avg = GHI_kt_avg;
		persistenceGHI.GHIavg = GHIavg;
		persistenceGHI.timeUTC = timeUTC;
		% Delete beginning entries from timeUTC, GHIavg, GHI_kt so that
		% everything lines up with the start time
		idf = find(timeUTC == starttime, 1);
		timeUTC(1:idf-1) = [];
		GHIavg(1:idf-1) = [];
		GHI_kt_avg(1:idf-1) = [];
	end
% 	% 30s averaging snippet
% 	for m = 1:floor(numel(GHIavg)/30)
% 		GHI30savg(m) = mean(GHIavg((m-1)*30+1:m*30));
% 		timeUTC30savg(m) = timeUTC(m*30);
% 	end
% 	% end snippet
else
	GHIavg=[]; GHI_kt_avg=[]; persistenceGHI=[]; GHIclrsky=[];
	warning('DEMROES:nositedata','There is no measured data for the selected site(s).');
end

end