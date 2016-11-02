function feeder = getFeeder(feederName, configId)
% this function is to generate different base configurations for the investigated feeder
%
% input:
%           feederName: could be name of the feeder or filename of the feeder config script.
%                       Excepted feeder config script names: Fallbrook.m or FallbrookConf.m
%           configId:   single or list of strings representing interested configs
%
% output:
%           feeders     :list of feeders in DSS struct format
%

% find corresponding feeder config file
feederName = strtrim(feederName);
if ~isempty(which(feederName))
    fncHandle = str2func(feederName);
elseif ~isempty(which([feederName 'Setup']))
    fncHandle = str2func([feederName 'Setup']);
elseif strcmp(feederName(end-1:end),'.m') % filename as input
    fncHandle = str2func(feederName);
else
    error('No config file for that feeder exist. Please double check the feederName/ path!');
end

feeder = fncHandle(configId);

end