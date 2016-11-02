function [fc, emptyProfileId] = fillForecastProfile(forecast,deploymentSite,fName,dayid)
% input is the forecast struct from the loadForecast function

emptyProfileId = [];
global conf; global indent;
if isempty(conf), conf = getConf; end
if ~exist('fName','var')
    d = getDeployment(deploymentSite);
    fName = d.Name; 
end
if ~exist('dayid','var'), dayid = ''; end
filepath = [conf.outputDir '/' fName '_Forecast_' dayid '_filled.mat'];
if exist(filepath,'file') 
    fprintf(['%sFilled forecast saved file exists! Load to use. File path: ' filepath '\n'],indent);
    fc = load(filepath); checkEmptyProfile(fc); return;
end

fc = forecast;
if ~exist('d','var'), d = getDeployment(deploymentSite); end

% footprint
fp = d.footprint;

% not nan/ empty data
avaiDat = ~isnan(fc.profile); 

% number of data points for each profile
ndat=sum(avaiDat);

% max ndat
ndatMax = max(ndat);

% fill each profile with neighbor's data until reaching ndatMax
for i = 1:length(ndat)
    if ndat(i) < ndatMax
        % get neighbor ids
        if isfield(fp,'pv')
            nb = findNearestSystem(fp.pv,i);
        else
            error('The findNearestSystem doesn''t handle footprint without pv struct in it yet. Write code to handle this. The code should use the footprint struct and pixel locations to figure out the nearest station');
        end
        
        nbid = 0;
        while ndat(i) < ndatMax
            nbid = nbid + 1;
            
            % fill them with the closest neighbor if available
            % data to use is the one that is empty/nan in this site and not nan/empty from the neighbor
            idToUse = isnan(fc.profile(:,i)) & ~isnan(fc.profile(:,nb(nbid)));
            
            % assign data
            fc.profile(idToUse,i) = fc.profile(idToUse,nb(nbid)); 
            
            % recalculate ndat(i)
            ndat(i) = sum(~isnan(fc.profile(:,i)));
        end
    end
end

save(filepath,'-struct','fc');
emptyProfileId = checkEmptyProfile(fc);
fprintf(['%sSaved filled forecast file: ' filepath '\n'],indent);
end

function d = getDeployment(deploymentSite)
% deployment site
if ischar(deploymentSite)
    d = siDeployment(deploymentSite);
else
    d = deploymentSite;
end
end