function fdata = loadforecastdata(imager, starttime, endtime, datatype, fields, datadir, ignoreSeriesFile)
% LOADFORECASTDATA load forecast results.
%
% data = loadforecastdata(imager,starttime,endtime,datatype,fields,datadir)
% loads the data of the given type with timestamps between starttime and
% endtime.
%
% Input:
%			starttime, endtime:
%						time bounds as datenum or a string datetime format.  Default endtime is now and default start time is 30 minutes prior to endtime
%			datatype:	one of the step type names.  Currently only 'forecast' has been tested
%			fields:		(optional) interested fieldname(s).  Default: all fields.  Note that many output files contain very large data and you will probably run out of memory if you try to load all the fields for a whole day at once.
%						Example: fields = {'time', 'ktavg'}
%			datadir:	(optional) forecast data directory in which each imager has its own folder named after its name.  Default is the outputDir param from forecast.conf
%

%% Process inputs
if ~exist('imager','var') || (~isa(imager,'imager') && ~ischar(imager))
	error('prettyplot:invalidInput','Please use a valid ''Imager'' object as input! Duhhh!');
elseif( ischar(imager) )
	imager = siImager(imager);
end
if ~exist('endtime','var') || isempty(endtime)
	endtime = nowUTC; % default to now
else
	if ischar(endtime)
		if( length(endtime) == 14 )
			endtime = datenum(endtime, 'yyyymmddHHMMSS');
		else
			endtime = datenum(endtime);
		end
	end
end
if ~exist('starttime','var') || isempty(starttime)
	starttime = endtime - 30/60/24; % default to 30 minutes from endtime
else
	if ischar(starttime)
		if( length(starttime) == 14 )
			starttime = datenum(starttime, 'yyyymmddHHMMSS');
		else
			starttime = datenum(starttime);
		end
	end
end
if ~exist('datatype','var') || isempty(datatype)
	datatype = 'forecast';
end
if ~exist('fields','var') || isempty(fields)
	fields = '*'; % load all
end
if ~exist('datadir','var') || isempty(datadir)
	conf = readConf(siGetConfPath('forecast.conf'));
	datadir = [conf.outputDir '/' imager.name];
else
    datadir = [datadir '/' imager.name];
end
if ~exist('ignoreSeriesFile','var') || isempty(ignoreSeriesFile)
	ignoreSeriesFile = 0;
end

%% get interested folder paths
DOYUTC = floor(starttime):floor(endtime);
inputDir = cell(length(DOYUTC),1);
for i = 1:length(DOYUTC)
	inputDir{i} = [datadir '/' datestr(DOYUTC(i),'yyyymmdd')];
end

%% look into each folder and find all the relevant files then concat them to a single variable
fdata = [];
for i = 1:length(inputDir)
	% load series file if exist
	if ~ignoreSeriesFile
		data = siLoadSeriesElmt(inputDir{i}, datatype, fields);
		if ~isempty(data)
			fdata = siJoinSeries(data,fdata);
			% refine time of interest
			id = (starttime < fdata.time(:,1)) & (fdata.time(:,1) < endtime);

			for fn = fieldnames(fdata)'; fn = fn{1};
				if iscell(fdata.(fn))
					fdata.(fn) = fdata.(fn){id};
				else
					fdata.(fn) = fdata.(fn)(id,:);
				end
			end
			if sum(id) > 0
				continue;
			end
		end
	end
	
	% load individual files 
	fl = dir([inputDir{i} '/' datatype '/' datatype '*.mat']);
	if ~isempty(fl)
		fl = {fl.name}';
		for j = 1:length(fl)
			fn = fl{j};
			% check if file is in the range of interest
			time = datenum(fn(length(datatype)+2:length(datatype)+15),'yyyymmddHHMMSS');
			if time >= starttime && time <= endtime
				if iscell(fields)
					data = load([inputDir{i} '/' datatype '/' fl{j}],fields{:});
				else
					data = load([inputDir{i} '/' datatype '/' fl{j}],fields);
				end
				fdata = siSeries(data,fdata);
			end
		end
	end
end

end