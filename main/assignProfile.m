function [c, prof] = assignProfile(c, type, prof, fcid, tzone, t, smoothFt)
% assign forecast profiles to targeted components (ES, PV, etc.)
% only supports pvsystem and load with daily profiles at the moment

global conf; global indent;
if isempty(conf), conf = getConf; end
% time for daily, need to change this somehow

if ~exist('t','var') || isempty(t)
    dt = conf.timeStep; t = 0 : dt/3600/24 : 1-dt/3600/24; t = t'; 
else
    dt = round( (t(2)-t(1))*24*3600 ); % in seconds.
end
% profile type
if isfield(conf,'mode'), profType = conf.mode; else profType = 'daily'; end

if ~exist('type','var') || isempty(type), type = 'pvsystem'; end
if ~exist('smoothFt','var'), smoothFt = 0; end

% convert nans in profiles to zeros
prof.profile(isnan(prof.profile)) = 0;


% % prof.profile(:,end+1:2*end)=prof.profile;
% % prof.profile(:,end+1:2*end)=prof.profile;

% convert profile timezone to simulation timezone if needed
if ~strcmpi(tzone,conf.timeZone)
    prof.time = timezoneConvert(prof.time, tzone, conf.timeZone);
end

% if the profile has been picked using fcid then set fcid to 1 to use that single profile for all systems
if length(fcid) == 1 && size(prof.profile,2) == 1
    fcid = 1;
end

% Theory: if fcid is defined, then use it as index for assigning
% elseif names are avail, use name of the targeted components for assigning data
% else assign assuming same indices are used in cir and fc
if ~exist('fcid','var') || isempty(fcid) || (length(fcid)==1 && ~fcid)
    % cut off some extra forecast profiles if not needed
    if size(prof.profile,2) >= length(c.(type)) % more forecast profiles than number of systems then use number of system as number of loadshapes
        prof.profile = prof.profile(:,1:length(c.(type)));
    else
        error('%snot enough forecast profiles!',indent);
    end
    lsid = 1:length(c.(type));
    
elseif length(fcid) == 1 || size(prof.profile,2) == 1% single profile for all systems
    prof.profile = prof.profile(:,fcid);
    lsid = ones(1,length(c.(type)));
    
elseif length(fcid) == length(c.(type)) % some specific ids are used
    prof.profile = prof.profile(:,fcid);
    lsid = 1:length(c.(type));
else
    error('%snumber of profiles are not the same as number of systems to be assigned to! Please check!',indent);
end

% convert profile with specified time to whole day profiles based on 't' (starting at midnight til 11:59:30 pm)
prof = makeWholeDayProfile(t,prof,dt,smoothFt,type);
lshape = createLoadShape(prof,dt,type);
c = addLShape(c,lshape);

% assign loadshape
for i = 1:length(c.(type))
    c.(type)(i).(profType) = lshape(lsid(i)).Name;
end

end

function lshape = createLoadShape(fc,dt,type)
lshape = dssloadshape;
for i = 1:size(fc.profile,2)
    lshape(i).Name = ['ls_' type  num2str(i)];
    lshape(i).Npts = length(fc.time);
    lshape(i).sInterval = dt;
    lshape(i).Mult = fc.profile(:,i);
end
end

function c = addLShape(c,lshape)
% add loadshapes to c
if isfield(c,'loadshape'), c.loadshape = [c.loadshape lshape];
else c.loadshape = lshape;
end
end