function [fc_, emptyProfileId] = loadForecast(dayid, fName, fcProfileId, forceReload)
% example of use: [fc,id] = loadForecast('20121219', 'fallbrook',45,1);
global conf; global indent;
if isempty(conf), conf = getConf; end
if ~exist('forceReload','var') || isempty(forceReload), forceReload = 0; end
fp = [conf.outputDir '/' fName '_Forecast_' dayid '.mat'];
if exist(fp,'file') && ~forceReload 
    fprintf(['%sForecast saved file exists! Load to use. File path: ' fp '\n'],indent);
    fc_ = load(fp); checkEmptyProfile(fc_); return;
end

% list forecast files 
fcfp = [conf.fcOutDir '/' conf.usi '/' dayid];
if exist(fcfp,'dir')
    p = dir([fcfp '/forecast_*']);
else
    fcfp = [conf.fcOutDir '/' fName '/' conf.usi '/' dayid];
    p = dir([conf.fcOutDir '/' fName '/' conf.usi '/' dayid '/forecast_*']);
end

if ~isempty(p)
p = sort({p.name}); p = p{end}; % sort by name, get latest updated file
fcfp = [fcfp '/' p];

	fc = load(fcfp);
	fnames = fieldnames(fc);
	
	switch conf.fcType
		case {'pv','inverter'}
			% find the right field name
			[v, id] = ismember({'pv','inverter'},lower(fnames));
			if any(v)
				if v(1) == 1
					fnId = 1;
				else
					fnId = 2;
				end
				type = fnames{id(fnId)};
			else
				error('%sThere is no ''pv'' or ''inverter'' field in the forecast output structure. Please double check!',indent);
			end
		otherwise
			error('%sHaven''t handle this case yet!',indent);
	end
	% clean up the forecast structure
	fieldsToRemove = setdiff(fnames,{type,'time'});
	fc = rmfield(fc,fieldsToRemove);
	
	% apply forecast horizon/ index of interest
	horId = conf.fcMin*2 + 1;
	fc_.time = fc.time(:,horId);
	if ~exist('opt','var') || ~isfield('opt','fcProfileId') || isempty(fcProfileId)
		for i = 1:size(fc.(type),2)
			x = [fc.(type)(:,i).gi];
			fc_.profile(:,i) = x(horId,:);
		end
	else
		% if specific profiles are wanted (saved as an array in opt.fcProfileId)
		for i = 1:length(fcProfileId)
			x = [fc.(type)(:,fcProfileId(i)).gi];
			fc_.profile(:,i) = x(horId,:);
		end
	end
	if isfield(conf,'GIfactor')
		GIfactor = conf.GIfactor;
	else
		GIfactor = 1000; % W/m2
	end
	fc_.profile = fc_.profile/GIfactor;
	save(fp,'-struct','fc_');
	fprintf(['%sSaved forecast file: ' fp '\n'],indent);
	emptyProfileId = checkEmptyProfile(fc_);
else
	fc_=[];
	emptyProfileId=[];
end