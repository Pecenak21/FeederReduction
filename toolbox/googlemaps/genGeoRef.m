function R = genGeoRef(coord)
	refCoord.lat(1) = min(coord.latitude(:));
	refCoord.lat(2) = max(coord.latitude(:));
	refCoord.lon(1) = min(coord.longitude(:));
	refCoord.lon(2) = max(coord.longitude(:));
	
	R = georasterref;
	R.RasterSize = size(coord.latitude);
	R.Latlim = refCoord.lat;
	R.Lonlim = refCoord.lon;
end