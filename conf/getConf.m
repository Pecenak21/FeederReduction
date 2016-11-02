function [conf, confPath] = getConf(varargin)
% obtain config parameters for simulation
%
% input:
%       simName     :simulation name & also folder name that contains the config info
%                   Location of this folder is specified in setup.conf (simSearchPath parameter)

%% get setup parameters
confPath = getConfPath('sim');

%% get config parameters for simulation from the config folder
conf = readConf( confPath ,1);
% handle string cell
for k = {'deployment', 'feederSetup','timeNote','timeDay','timeStart','timeEnd','loadProfTime','feederName'}
    i = k{1};
    conf.(i) = strtrim(strsplit(conf.(i),','));
end

% handle file path
for k = fieldnames(conf)'
    i = k{1};
    if ~iscell(conf.(i)) && ~isempty(strfind(conf.(i),'$'))
        conf.(i) = normalizePath(conf.(i));
    end
end


% handle num cell
for k = {'PVLevel'}
    i = k{1};
    if ischar(conf.(i))
        conf.(i) = str2double(strtrim(strsplit(conf.(i),',')));
    end
end

% handle time
% conf.timeStep = conf.timeStep/24/3600; % convert to day
% for k = {'timeStart','timeEnd','timeDay','loadProfTime'}
%     i = k{1};
%     conf.(i) = datenum(strtrim(strsplit(conf.(i),',')));
%     % convert to UTC if they are not in UTC
%     if ~strcmpi(strtrim(conf.timeZone),'utc')
%         try conf.(i) = timezoneConvert(conf.(i),strtrim(conf.timeZone),'UTC');
%         catch e
%             error(['Not supported timezone: ' conf.timeZone '. Avail time zones can be found by running java.util.TimeZone.getAvailableIDs in Matlab or in the timeZoneList.m file in gridIntegration repo.']);
%         end
%     end
% end

%% gridIntegration root dir
conf.gridIntDir = pwd;

%% add other parameters if given
if length(varargin)>1 && mod(nargin,2) == 0
    vars = reshape(varargin, 2, nargin/2);
    for i = 1:size(vars,2)
        try
            conf.(vars{1,i}) = vars{2,i};
        catch
            error('Field name must be a string!');
        end
    end
end

end