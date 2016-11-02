function [d, dat] = getLoadProfile(fName,loadProfileId,tid)
% Load the data and select the profiles of interest
% currently only support daily profile when specified
%
% output:
%           d   : filtered data
%           dat : raw data

global conf; global indent;
if isempty(conf), conf = getConf; end

fp = [conf.outputDir '/' fName '_loadProfile.mat'];
% if exist(fp,'file')
%     fprintf(['%sLoad profile''s saved file exists! Load to use. File path: ' fp '\n'],indent);
%     x = load(fp); d = x.d; dat = x.dat; return;
% end

% load the file
d = load(conf.loadProfile); dat = d;

% load profile of interest
if loadProfileId
    d.profile = d.profile(:,loadProfileId);
    d.profileNames = d.profileNames(loadProfileId);
end

if isfield(conf,'loadProfTime') && ~isempty(conf.loadProfTime)
    conf2 = conf;
    conf2.loadProfTime=conf.timeDay(tid);
    d = getProfile(d, conf2);
end

if conf.loadProfScaling
    lmin = 0; lmax = 1;
    if conf.loadProfMin, lmin = conf.loadProfMin; end
    if conf.loadProfMax, lmax = conf.loadProfMax; end
    dmin = repmat(min(d.profile),size(d.profile,1),1);
    dscaling = repmat((lmax - lmin)./ (max(d.profile) - min(d.profile)),size(d.profile,1),1);
    
    d.profile = (d.profile - dmin) .* dscaling + lmin;
end

save(fp,'d','dat');
end

function d = getProfile(d, conf)
time = d.time; prof = d.profile; d.time = []; d.profile = [];
% convert load's timezone if needed to local time for simulation
if ~strcmpi(conf.timeZone,conf.loadProfTZone)
    time = timezoneConvert(time,conf.loadProfTZone,conf.timeZone);
end

% get relevant time index out of the data series based on the time periods/days specified
tid = [];
for i = 1:length(conf.loadProfTime)
    t = datenum(conf.loadProfTime{i});
    tid = [tid find((time >= t) & (time < t+1 ))];
end

% get data
d.time(:,1) = time(tid);
d.profile(:,:) = prof(tid,:);

end