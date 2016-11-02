%%
orig_coord = '37.869733350860730,-122.284758688533';
dest_coord = '37.871700000000000,-122.253300000000';
mode='walking';

url = ['https://maps.googleapis.com/maps/api/distancematrix/json?origins=(',orig_coord,')&destinations=(',dest_coord,')&mode=',mode,'&language=en-EN&sensor=false'];

str = urlread(url);

https://maps.googleapis.com/maps/api/distancematrix/json?origins=(37.869733350860730,-122.284758688533)&destinations=(37.871700000000000,-122.253300000000)&mode=walking&language=en-EN&sensor=false
https://maps.googleapis.com/maps/geo?output=xml&key=AIzaSyD1tUkziDpDsb8IkmyVZmHbydBjyp-G174&q=ucsd

%%
filenames = gunzip('sanfranciscos.dem.gz', tempdir); 
demFilename = filenames{1}; 

[lat, lon,Z] = usgs24kdem(demFilename,2); 

% Delete the temporary gunzipped file. 
delete(demFilename); 

% Move all points at sea level to -1 to color them blue. 
Z(Z==0) = -1;

% Compute the latitude and longitude limits for the DEM. 
latlim = [min(lat(:)) max(lat(:))];
lonlim = [min(lon(:)) max(lon(:))];

% Display the DEM values as a texture map. 
figure
usamap(latlim, lonlim)
geoshow(lat, lon, Z, 'DisplayType','texturemap')
demcmap(Z)
daspectm('m',1)

% Overlay black contour lines onto the texturemap.
geoshow(lat, lon, Z, 'DisplayType', 'contour', ...
  'LineColor', 'black');

%%
   lat = [48.8708   51.5188   41.9260   40.4312   52.523   37.982];
   lon = [2.4131    -0.1300    12.4951   -3.6788    13.415   23.715];
   plot(lon,lat,'.r','MarkerSize',20)
   plot_google_map