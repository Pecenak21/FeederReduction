function tl = LinkeTurbidity( time, latitude, longitude )
% handle input
if(isstruct(latitude)) % user gave a position struct instead of lat/lon/alt
	longitude = latitude.longitude;
	latitude = latitude.latitude;
end
% persistent cache of linke data
persistent linkeDB;
% Linke Turbidity database is stored by month
month = datevec(time); month = month(:,2);
[months, ~, im] = unique(month);

% fill in any unloaded fields in our db
if(isempty(linkeDB)), linkeDB = cell(1,12); end
if(any(cellfun('isempty',linkeDB(months))))
	conf = readConf(getConfPath('setup.conf'));
	d = normalizePath(conf.LINKE_PATH);
	for m = months(:)'
		if(isempty(linkeDB{m}))
			linkeDB{m} = linkeRead([d '/TL5_' lower(datestr([2000 m 1 0 0 0], 'mmm')) '.bin']);
		end
	end
end
% scale the latitude and longitude appropriately: 1degree = 60' = 12 pixels
latitude = round((-latitude+90)*12 + 0.5);
longitude = round((longitude+180)*12 + 0.5);

% These three lines would do it the way Bryan's java code does it, which I claim is probably not right...
% ss = 360/(4320-1);
% latitude = floor((90-latitude)/ss+1+1);
% longitude = floor((longitude+180)/ss+1+1);

s = size(linkeDB{months(1)});
% this scaling puts 90S and 180E on a pixel off the edge of the map, so we need to bring them back on
% note that this is reasonable because (90-1e9)S
if(latitude > s(1)), latitude = s(1); end
if(longitude > s(2)), longitude = s(2); end
% look up the turbidity for the desired months at the given location
tl = cellfun(@(x)x(latitude,longitude),linkeDB(months));
% copy across all times
tl = reshape(tl(im),size(time));
end

function dat = linkeRead(file)
% read raw binary data from the linke db file.  From the readme:
%
% Gridded data, raw format, 2160 rows and 4320 columns, 1 byte per value. Cell size is 5'.
% Upper left corner is 90 N, 180 W. Then, point 90 N, 179.5 W etc. Lower right is 90S 180 E.
%
% The Linke turbidity factor has no unit. It typically ranges between 3
% (clear skies) to 7 (heavily polluted skies). The factor was multiplied by
% 20 for storage. A value of 70 means a Linke turbidity factor equal to
% 3.5.

fid = fopen(file);
dat = fread(fid, [4320, 2160])'/20; %matlab wants to fill down columns rather than along rows
fclose(fid);

end
