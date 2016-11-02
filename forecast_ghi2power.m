function [fc_, emptyProfileId] = forecast_ghi2power(dayid, fName, fcProfileId, forceReload,Input_data)
%This function is a copy of loadforecast.m with minor modifications to
%adjust loads more accurately using the Perf model. The Additional input
%variable "Input_data" is an optional structure that can have fields that
%may be passed to the Perf model to either change default values used in
%the Perf code, or to do additional adjustments, including adjustments temperature and
%wind
%
%In this funcation the main change of the code is the loading of GHI
%instead of gi and the calling of the perf_model


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
p = sort({p.name}); p = p{end};%p = 'forecast_002.mat'; % sort by name, get latest updated file
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
        x = [fc.(type)(:,i).ghi];
        y = [fc.(type)(:,i).ktavg];
        fc_.profile(:,i) = x(horId,:);
        fc_.kt(:,i) = y(horId,:);
    end
else
    % if specific profiles are wanted (saved as an array in opt.fcProfileId)
    for i = 1:length(fcProfileId)
        x = [fc.(type)(:,fcProfileId(i)).ghi];
        y = [fc.(type)(:,fcProfileId(i)).ktavg];
        fc_.profile(:,i) = x(horId,:);
        fc_.kt(:,i) = y(horId,:);
    end
end

%% feeder id
[~,fid] = ismember(lower(fName),lower(conf.feederName));

%% Create required input for Perf model and call it
target=siDeployment(conf.deployment{fid});
Input_data.time=fc_.time;
Input_data.GHI=fc_.profile;
Input_data.kt=fc_.kt;
time_diff=datenum([0 0 0 23-conf.fcTimeAhead 0 0]);
[temp, wind]=load_interp_QCLCD(Input_data.time-time_diff,conf.feederName{fid});

if ~isempty(temp)
    Input_data.Temp_Amb=temp;
end

if ~isempty(wind)
    Input_data.Wind_Speed=wind;
end
[power]=Perf_model(target,Input_data);

%Place power data in output and save
fc_.profile=power;
if isfield(conf,'GIfactor')
    GIfactor = conf.GIfactor;
else
    GIfactor = 1000; % W/m2
end
fc_.profile = fc_.profile/GIfactor;
save(fp,'-struct','fc_');
fprintf(['%sSaved forecast file: ' fp '\n'],indent);
emptyProfileId = checkEmptyProfile(fc_);
end