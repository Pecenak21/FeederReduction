classdef pvsystem < handle
	%PVSYSTEM Summary of this class goes here
	%   Detailed explanation goes here
	
	properties (SetObservable, GetObservable)
		name;
		cenlat; % center latitude
		cenlon; % center longitude
		lat; % lat = [minlat maxlat]
		lon; % lon = [minlin maxlon]
		area;
		kVA;
		kW;
		dlat; % 
		dlon; % 
		dxm; % width in meter
		dym; % height in meter
		outputrating = 100; % W/m^2
		zenTilt = 0; % in degrees. Default: flat
		azTilt = 180; % in degrees. Default: South oriented
		efficiency = .1; % 10% 
		density = .95; % how packed the PV panels are arranged
		footprint; % map with lat lon limit based on lat lon params
		note;
	end
	
	methods
		function p = calArea( p )
			% calculate the area of PV based on its footprint defined by
			% lat and lon params. If lat and lon are not known, estimate
			% the the required area based on PV location (cenlat,cenlon),
			% PVsize in kVA, and density of PV system.
			
			if isempty(p.lat) || isempty(p.lon)
				p.calAreaVirtual;
				return;
			end
			
			if isempty(p.cenlat) || isempty(p.cenlon)
				p.cenlat = mean(p.lat);
				p.cenlon = mean(p.lon); 
			end
			
			p.calAreaVirtual(1);
			
			p.dym = 1000*greatCircleDistance( p.lat(1), mean(p.lon),...
									p.lat(2), mean(p.lon) );
			p.dxm = 1000*greatCircleDistance( mean(p.lat), p.lon(1),...
									mean(p.lat), p.lon(2) );
			
			% calculate density
			p.density = p.area/(p.dym*p.dxm);
			
			% area
			p.area = p.dxm * p.dym;
			
			% dlat, dlon
			p.dlat = p.lat(2) - p.lat(1);
			p.dlon = p.lon(2) - p.lon(1);
		end
		
		%% organize lat data when it is set. 
		function set.lat(obj,val)
			obj.lat = [min(val(:)) max(val(:))];
		end
		
		%% organize lon data when it is set. 
		function set.lon(obj,val)
			obj.lon = [min(val(:)) max(val(:))];
		end
	end
	
	methods (Access=private)
		function p = calAreaVirtual( p, density )
			% If defined footprint (from lat, lon params) is not known,
			% estimate the area based on PV size and center location (cenlon,cenlat). 
			
			if exist('density','var') && ~isempty(density)
				p.density = density;
			end
			
			% area needed [m^2] = size [kW] * 1000 / density [1] / outputrating [W/m^2]
			p.area = p.kVA * 1000 / p.density / p.outputrating;

			% size of the square
			p.dxm = sqrt(p.area); p.dym = p.dxm;

			%% how much is a distance (in m) of .01 change in latitude at PV location
			dl = .01;
			d = 1000*greatCircleDistance(p.cenlat - dl/2, p.cenlon,...
									p.cenlat + dl/2, p.cenlon);
			% dlat needed
			p.dlat = dl/d*p.dxm;

			% dlon needed
			d = 1000*greatCircleDistance(p.cenlat, p.cenlon - dl/2,...
									p.cenlat, p.cenlon + dl/2);
			p.dlon = dl/d*p.dxm;
			
			% lat, lon
			if isempty(p.lat) || isempty(p.lon)
				p.lat = [p.cenlat - p.dlat/2; p.cenlat + p.dlat/2];
				p.lon = [p.cenlon - p.dlon/2; p.cenlon + p.dlon/2];
			end
		end
	end
	
end

