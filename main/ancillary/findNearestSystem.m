function [nearestSystemId, d] = findNearestSystem(dat,systemId,excludeSystem)
% dat is a struct that has 'cenlat' and 'cenlon' fields for central latitude and longitude position

lat = dat(systemId).cenlat;
lon = dat(systemId).cenlon;

% find distance from system of interest to every other system
d = sqrt((lat - [dat.cenlat]).^2 + (lon - [dat.cenlon]).^2);

% assign the distance to self and excluding systems to inf
d(systemId) = inf;
if exist('excludeSystem','var') && ~isempty(excludeSystem)
    d(excludeSystem) = inf;
end

[~,nearestSystemId] = sort(d);
end