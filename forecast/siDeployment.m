function d = siDeployment(dName)
% siDeployment loads all the configuration options needed to forecast for a deployment site
% The returned deployment is a struct with appropriate elements.
%
% siDeployment(name) looks for a deployment with the given name in the search path specified in siDeployment.conf
% siDeployment(dir_path) loads the deployment info from the specified directory
%
% A few earlier deployments may not conform to this structure, but the loaded structure should contain the following:
%	d.name - filled from folder name
%	d.data_type - cell array of strings (or single string) giving the names of footprints  present.  Typically 'PV', or 'sensor', though others could conceivably be added.
%	d.data_fcn - cell array of strings (or single string) describing functions that get timeseries data for each data type, for pv/ghi data from different sources
%		* data_fcn will be called once for each data type, as
%		  [rawdata, rawtime] = data_fcn(start_time, end_time, target);
%		  rawdata and rawtime should be a cell array of column vectors
%		  containing data values and the corresponding times for each
%		  sensor.  It is suggested that you simply return all available
%		  data in the time window and allow the forecast code to
%		  interpolate as desired.  The cell array format allows each
%		  sensor/pv system to have 'available data' at different times.
%	d.ground - structure containing rectangular grid of latitude, longitude, and altitude of the ground at the site
%	d.footprint - structure containing image maps that show which pixel in the groundmap correspond to each sensor/PV array
%				  footprint also generally contains names and input data used to generate the respective maps
%	d.design - structure containing info about tilt (positive, away from zenith), azimuth (degrees clockwise from North), and nominal power (W) of the item of interest
%		The following special considerations apply to these fields:
%		* specify tilt,azimuth = -1 for 2-axis tracking (e.g. DNI sensor).  single axis tracking is not currently defined.  Not currently implemented in forecast
%		* specify multiple values in each row for nominal power to have the values passed to polyval instead of multplied directly (polynomial coefficients should give power in watts from modeled GI in kW)
%		* for GHI/DNI, nominalpower should be 1000 to convert from kW to W.
%
% See Also: siDeploymentTool, siForecastGHI, siPowerElmt, str2func

%% Start by loading the main conf file for this deployment
depname = dName;
if( ~exist([dName '/deployment.conf'],'file') )
	searchPath = readConf(getConfPath('setup.conf'));
	searchPath = searchPath.siDeploymentSearchPath;
	if( ischar(searchPath) )
		searchPath = siNormalizePath(searchPath);
		if( exist([searchPath '/' dName],'dir') )
			dName = [searchPath '/' dName];
		end
	elseif iscell(searchPath)
		for s = searchPath(:)'
			s = normalizePath(s{1});
			if( exist([s '/' dName],'dir') )
				dName = [s '/' dName];
				break;
			end
		end
	end
	if( ~exist(dName,'dir') )
		error('Could not find a deployment named "%s"',dName);
  end
else
  dName = normalizePath(dName);
end
%% Parse the conf file
d = readConf([dName '/deployment.conf']);

d.name = depname;
flist = [{'name'};fieldnames(d)];
d = orderfields(d,flist(1:end-1));
d.source = dName;

% load external files
flist = dir([dName '/*.mat']);
for f = flist'
  try
    d.(regexprep(f.name,'\.mat$','')) = load([dName '/' f.name]);
  catch e %#ok<NASGU>
    warning( ['File ' dName '/' f.name ' did not load.' ] ); %#ok<WNTAG>
  end
end

%% Specialized parsing of deployment file


% Convert latitude, longitude and altitude to a single substruct
if( sum(isfield(d,{'longitude';'latitude';'altitude'})) == 3 )
	d.position.longitude = d.longitude; d = rmfield( d , 'longitude' );
	d.position.latitude  = d.latitude;  d = rmfield( d , 'latitude' );
	d.position.altitude  = d.altitude;  d = rmfield( d , 'altitude' );
end

% If we already have a groundmap, delete the parameters that would be used to generate it just to keep things tidy
if( isfield(d,'ground') )
	fn = {'position';'gmap_n';'gmap_e';'gmap_s';'gmap_w';'gmap_res';'gmap_generate';'gmap_save'};
	fn( ~isfield( d , fn ) ) = [];
	d = rmfield(d,fn);
end

% if we don't have a groundmap, go ahead and generate one and save it
if( ~isfield(d,'ground') && all(isfield(d,{'gmap_n';'gmap_e';'gmap_s';'gmap_w';'gmap_res';'position'})) )
	% generate a groundmap of the desired size
	nesw = [ d.gmap_n d.gmap_e d.gmap_s d.gmap_w ];
	d.ground = geo_getGroundmap( d.position , nesw , d.gmap_res );
	% remove the extra fields that we don't need
	d = rmfield(d,intersect({'gmap_n';'gmap_e';'gmap_s';'gmap_w';'gmap_res';'gmap_generate';'gmap_save';'position'},fieldnames(d)));
	% save for use the next time
	ground = d.ground; %#ok<NASGU>
	save( [dName '/ground.mat'] , '-struct' , 'ground' );
end

% Load any pi data sources
if( isfield(d,'PI_CONF') )
	d.pi.obj = piServer(siGetConfPath(d.PI_CONF));
end

if( isfield(d,'UCSDPI_CONF') )
	d.ucsdpi = ucsdPiServer;
end

%% legacy data type support
% Henderson-era deployments may be missing data_type field.  Value is 'inverter'
if(~isfield(d,'data_type'))
	d.data_type = 'inverter';
end
% UCSD augmented deployment 
if(~iscell(d.data_type) && strcmpi(d.data_type,'ghipv'))
	d.data_type = {'GHI','pv'};
	d.footprint.pvnames = d.footprint.inverterNames;
	d_pv = d.footprint.pv; % preserve PV design specs
	d.footprint.pv = d.footprint.GHI;
	d.footprint.pv(d.footprint.pv < 100) = nan;
	d.footprint.pv = d.footprint.pv - 100;
	d.footprint.GHI(d.footprint.GHI>100) = nan;
	% update design structure
	d.design.pvtilt = d_pv.tilt;
	d.design.pvazimuth = d_pv.azimuth;
	% old form was an efficiency coefficient, new form is applied to power directly
	p_old = [d_pv.eff_coeff 0];
	d.design.pvnominal = p_old .* ((d_pv.scaleFactor*1000).^((numel(p_old)-1):-1:0));
end
% make sure the data type is a cell now
if(~iscell(d.data_type)), d.data_type = {d.data_type}; end

% add data functions to UCSD deployments, using old condition from siPowerElmt
if ~isempty(regexpi(d.name,'^ucsd')) && isfield(d,'data_type') && strcmpi(d.data_type{1}(1),'g') && ~isfield(d,'data_fcn')
	d.data_fcn = {@getDEMROES};
	if(numel(d.data_type)>1)
		d.data_fcn{2} = @(st,et,t)t.ucsdpi.getPower(st,et,1);
	end
% elseif(isfield(d,'data_fcn'))
% 	if(ischar(d.data_fcn))
% 		d.data_fcn = str2func(d.data_fcn);
% 	elseif(iscellstr(d.data_fcn))
% 		d.data_fcn = cellfun(@str2func,d.data_fcn, 'UniformOutput', false);
% 	end
end

%% if any power/sensor data sources have polynomial calibrations, load the inverse functions
for di = 1:numel(d.data_type);
	dt = d.data_type{di};
	% add a default nominal field if it doesn't exist
    try 
        s = [numel(d.footprint.([dt 'names'])),1];
    catch
        s = [numel(d.footprint.([dt 'Names'])),1];
    end
	if(~isfield(d.design,[dt 'nominal']))
		d.design.([dt 'nominal']) = 1000*ones(s);
	end
	if(~isfield(d.design,[dt 'tilt']))
		d.design.([dt 'tilt']) = zeros(s);
		d.design.([dt 'azimuth']) = zeros(s);
	else
		if(~iscolumn(d.design.([dt 'tilt'])))
			warning('siDeployment:tilt_as_column','tilt and azimuth should be specified as column vectors.  If they aren''t many things may break.  Here, let me fix that for you (but you should do it in the file to avoid future warnings).');
		end
	end
	m = any(d.design.([dt 'nominal'])(:,1:end-1),2);
	d.design.([dt 'nominal_is_poly']) = m;
	if(any(m))
		x = (-.05:.01:1.6)'; % valid range in kW/m^2 for power calibration
		d.design.([dt 'nominal_inv']) = cell(s(1),1);
		coeff = d.design.([dt 'nominal']);
		for fi=find(m)
			y = polyval(coeff(fi,:),x);
			% may need to reduce the range slightly depending on the given polynomial
			ymask = [false; sign(diff(y))>0]; ymask(diff(ymask)>0) = true;
			if(~issorted(y(ymask)))
				warning('siDeployment:calibrationPoly','bad calibration polynomial may be too far from monotonic for me to correct');
				ydown = find(diff(ymask)<0);
				yup = find(diff(ymask)>0);
				[~,imax] = max(ydown-yup);
				ymask(1:yup(imax)) = false; ymask((ydown(imax)+1):end) = false;
			elseif(min(y(ymask))>0 || max(y(ymask))<1.3)
				warning('siDeployment:calibrationPoly','Calibration polynomials should be monotonically increasing at least from 0 to 1.3.');
			end
			d.design.([dt 'nominal_inv']){fi} = griddedInterpolant(y(ymask),x(ymask));
		end
	end
end

end
