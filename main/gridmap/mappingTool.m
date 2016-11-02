function h = mappingTool(powerGridMapPath)
% powerLineMapping maps 1 line diagram to geo-coordinate map (using google
% map)

%% 1. load power grid map
%% draw figure window
f = figure;
ss = get(0,'screensize');
set(f,'toolbar','figure','position',[50, 10, ss(3:4)-100]);


%% Plot power grid map
% Map
h.axes = axes('position',[0.03 0.03 0.77 0.94]);
h.map = imshow(powerGridMapPath); %plot_google_map('MapType', 'hybrid');

%% control panel
% Deployment Name
uicontrol('style','text','units','normalized', 'position', [0.82, 0.95, 0.16, 0.02], 'String', 'Deployment Name:');
h.name = uicontrol('style','edit', 'units','normalized', 'position', [.82, .92, .16, .03], 'String', 'Default Deployment Name');
% Resolution
uicontrol('style','text','units','normalized','position', [.82, .87, .10, .04], 'String', 'Forecast Resolution (meters/pixel)');
h.resText = uicontrol('style','edit','units','normalized','position', [.92, .87, .06, .04], 'String', '2.5');
h.resSlider = uicontrol('style','slider', 'units', 'normalized', 'position', [.82, .85, .16, .02], 'Min', .1, 'Max', 10, 'Value', str2double(get(h.resText,'string')), 'Callback', @(obj,ev)set(h.resText,'string',num2str(get(obj,'value'))));
set(h.resText, 'callback', @(obj,ev)set(h.resSlider,'value',str2double(get(obj,'string'))));
% Extra space for border
uicontrol('style','text','units','normalized','position', [.82, .81, .10, .03], 'String', 'Border size (%)');
h.borderText = uicontrol('style','edit','units','normalized','position', [.92, .81, .06, .03], 'String', '10');
% Table of buses
uicontrol('style','text','units','normalized', 'position', [0.82, 0.78, 0.16, 0.02], 'String', 'Buses:');
h.buses = uitable('Units', 'normalized', 'position', [.82 .48 .16 .3], 'ColumnName', {'Name','Lat','Lon'}, 'ColumnFormat', {'char','numeric','numeric'}, 'columnEditable', true, 'columnWidth',{100,60,60}, 'CellSelectionCallback', @tableCellSelect);
h.busAdd = uicontrol('units','normalized','position', [.9, .46, .02, .02], 'string', '+', 'Callback', @addBus);
h.busDel = uicontrol('units','normalized','position', [.92, .46, .02, .02], 'string', '-', 'Callback', @rmBus);
h.calibrate = uicontrol('units','normalized','position', [.94, .46, .04, .02], 'string', 'Calibrate', 'Callback', @calibrate);
% Table of lines
uicontrol('style','text','units','normalized', 'position', [0.82, 0.43, 0.16, 0.02], 'String', 'Lines:','enable','off');
h.lines = uitable('Units', 'normalized', 'position', [.82 .13 .16 .3], 'ColumnName', {'Name','Bus1','Bus2'}, 'ColumnFormat', {'char','numeric','numeric'}, 'columnEditable', true, 'columnWidth',{120,50,50}, 'CellSelectionCallback', @tableCellSelect);
set(h.lines,'enable','off');
h.lineAdd = uicontrol('units','normalized','position', [.94, .11, .02, .02], 'string', '+', 'Callback', @addLine);
set(h.lineAdd,'enable','off');
h.lineDel = uicontrol('units','normalized','position', [.96, .11, .02, .02], 'string', '-', 'Callback', @rmLine);
set(h.lineDel,'enable','off');
% Save buttons
h.save = uicontrol('units', 'normalized', 'position', [.82, .02, .075, .05], 'string', 'Save to File', 'Callback', @saveDeployment);
h.saveToWorkspace = uicontrol('units', 'normalized', 'position', [.905, .02, .075, .05], 'string', 'Save to Workspace', 'Callback', @saveDeployment);

% save the control handles with the figure data
setappdata(f,'Controls',h);

%% 2. pick first point in the map and give lat lon info (obtained from google map for example)
msgbox('Please pick at least 2 points/ buses and input their lat-lon coordinate info for calibration');
global xfit; xfit = [];
global yfit; yfit = [];
% if(nargin > 0)
% 	% draw point of imager on map
% 	hold on; plot(position.longitude, position.latitude', 'rx'); hold off;
% end
%% 3. pick 2nd point in the map and give lat lon info (obtained from google map for example)



%% 6. when ready to go, ask to pick the first bus by left click
% after picking, pop up a window asking for name of the bus, if not, use
% 'b1'
%% 7. after that, for all subsequent buses, as if there is a line between previous bus and the one just added

%% 8. allow user to choose adding line between established buses (do add bus and add line), allow clicking 2 buses to make line

%% 9. allow user to add line properties

%% 10. save line and bus info

%% calibrate
function calibrate(obj, ev)
button = questdlg('Have you input the lat-lon info for the known buses?','Reminder for lat-lon info before calibration','No','Yes','Yes');
switch button
    case 'No'
        msgbox('Calibration cancelled!'); return;
end
f = get(obj,'parent');
h = getappdata(f,'Controls');
d = get(h.buses,'Data');

%% overlay power grid map on google map to see if the input is correct
% calculate the boundaries for google map 
p = getappdata(f,'BusCoordinates');

xfit = polyfit(p(:,1)',[d{:,2}],1);
yfit = polyfit(p(:,2)',[d{:,3}],1);

% xbound = sort(polyval(xfit,get(h.map,'XData')));
% ybound = sort(polyval(yfit,get(h.map,'YData')));
% 
% d = get(h.axes,'position')
% 
% axis([ybound, xbound]);
% 
% h.map = plot_google_map('MapType', 'roadmap');
% set(h.map,'alphadata',0.6);
% hold on
% im = imread(powerGridMapPath);
% imshow(powerGridMapPath);
% h.map2 = imshow(powerGridMapPath);
% plot 
% position = struct('latitude',0,'longitude',0,'altitude',0);
% if(nargin>0)
% 	if(ischar(imager))
% 		imager = siImager(imager);
% 	end
% 	position = imager.position;
% else
% 	LatBounds = [-70, 70];
% 	LonBounds = [-180, 180];
% end
% get the size of a degree at the central position (need to update this later again)
% [meters_degLat, meters_degLon] = sizeOfADegree(position.latitude);
% % assign bounds for if we have a central imager; assume 5km radius for now
% if(~exist('LatBounds','var'))
% 	LatBounds = position.latitude + [-1,1]*5000/meters_degLat;
% 	LonBounds = position.longitude + [-1,1]*5000/meters_degLon;
% end
% axis([LonBounds, LatBounds]);
%% 5. ask if the input is correct, if not (want to modify) go back to step 2

% disp('Does it look good?'); 
% if yes, then done and disable after calibrating
%set(h.calibrate,'enable','off');

msgbox('Calibration done!');

end

end

% Functions

% Add new polygon
function addLine(obj, ev)
f = get(obj,'parent');
h = getappdata(f,'Controls');

% have the user draw the new polygon
figure(f);
title('Trace line.  When done, right click and select "Create Line"');
[~, xi, yi] = roipoly();
title('');

if(isempty(xi)), return; end

% save the coordinates of the new polygon
p = getappdata(f,'PVcoordinates');
if(isempty(p))
	p = {[xi,yi]};
else
	p(end+1,:) = {[xi,yi]};
end
setappdata(f,'PVcoordinates',p);

% draw the new polygon
hh = patch(xi,yi,size(p,1));
set(hh,'Tag', sprintf('pv%03d',size(p,1)));

% add the new region to the list:
% for convenience, we try to guess names of the form 'xxx 1', 'xxx 2', etc.
names = get(h.lines, 'Data');
if(isempty(names))
	root = [];
else
	root = regexp(names{end,1},'(.*?)(\d+)$','tokens','once'); % ? makes it lazy instead of greedy
end
if(~isempty(root))
	names{end+1,1} = sprintf('%s%i',root{1}, str2double(root{2})+1);
else
	names{end+1,1} = sprintf('Line %i', size(p,1));
end
names(end,2:4) = {0,180,1000};
set(h.lines, 'Data', names);

end

% Add new point
function addBus(obj, ev)
f = get(obj,'parent');
h = getappdata(f,'Controls');

% Have the user input a new point
figure(f);
title('click on location of the bus');
set(f, 'WindowButtonDownFcn', @(obj,ev)finishAddBus(obj,ev,f,h));
% waitfor(h.axes,'CurrentPoint');
end

function finishAddBus(~,~,f,h)
global xfit;
global yfit;
set(f, 'WindowButtonDownFcn', []);
x = get(h.axes,'CurrentPoint');
y = x(1,2); x = x(1,1);

if(isempty(x)), return; end

% save the new coordinates
p = getappdata(f,'BusCoordinates');
if(isempty(p))
	p = [x,y];
else
	p(end+1,:) = [x,y];
end
setappdata(f,'BusCoordinates',p);

% draw the point
figure(f); hold on;
hh = plot(x,y,'s', 'linewidth',2);
hold off;
set(hh,'Tag', sprintf('pt%03d',size(p,1)));
title('');

% add the new region to the list
d = get(h.buses, 'Data');
if(isempty(d))
	root = [];
else
	root = regexp(d{end,1},'(.*?)(\d+)$','tokens','once'); % ? makes it lazy instead of greedy
end
if(~isempty(root))
	d{end+1,1} = sprintf('%s%i',root{1}, str2double(root{2})+1);
else
	d{end+1,1} = sprintf('Bus %i', size(p,1));
end
d{end,2} = x; d{end,3} = y;
if ~isempty(xfit)
    d{end,2} = polyval(xfit,d{end,2});
    d{end,3} = polyval(yfit,d{end,3});
end
set(h.buses, 'Data', d);

end

% remove point or polygon
function rmLine(obj, ev)
f = get(obj,'parent');
h = getappdata(f,'Controls');
% delete the item from the table
ind = getappdata(h.lines, 'selectedcell');
if(isempty(ind)), error('no polygon selected'); end
ind = ind(1);
names = get(h.lines, 'Data');
names(ind,:) = [];
set(h.lines, 'Data', names);
% delete the patch for this item
patches = findobj(h.axes, 'type', 'patch', 'tag', sprintf('pv%03d', ind));
delete(patches);
% decrement the "color" of all subsequent patches by 1
patches = findobj(h.axes, 'type', 'patch');
for i = 1:numel(patches);
	x = get(patches(i),'tag');
	if(numel(x)<2 || ~strcmp(x(1:2),'pv')), continue; end
	x = str2double(x(3:end));
	if(x > ind)
		set(patches(i), 'cdata', x-i);
		set(patches(i), 'tag', sprintf('pv%03d', x-i));
	end
end
% remove the vertex data from that array
polygons = getappdata(f,'PVcoordinates');
polygons(ind) = [];
setappdata(f,'PVcoordinates',polygons);
end

function rmBus(obj, ev)
f = get(obj,'parent');
h = getappdata(f,'Controls');
% delete the item from the table
ind = getappdata(h.buses, 'selectedcell');
if(isempty(ind)), error('no point selected'); end
ind = ind(1);
names = get(h.buses, 'Data');
names(ind,:) = [];
set(h.buses, 'Data', names);
% delete point for this item
points = findobj(h.axes, 'type', 'line', 'tag', sprintf('pt%03d', ind));
delete(points);
% decrement the tag for all subsequent points by 1
points = findobj(h.axes, 'type', 'line');
[mask, pointnums] = regexp(get(points,'tag'), 'pt(\d{3})','once','start','tokens');
if(iscell(mask))
	mask = cellfun(@isempty,mask);
else
	mask = isempty(mask);
end
points(mask) = [];
pointnums = str2double([pointnums{~mask}]);
mask = pointnums<ind;
points(mask) = [];
pointnums(mask) = [];
for i = 1:numel(points)
	set(points(i),'tag', sprintf('pt%03d', pointnums(i)-1));
end
% remove the point data:
p = getappdata(f, 'BusCoordinates');
p(ind,:) = [];
setappdata(f, 'BusCoordinates', p);
end

% keep track of selected index on tables
function tableCellSelect(obj,ev)
setappdata(obj, 'selectedcell', ev.Indices);
end

% Rename point or polygon
% I think this can be implemented just letting the table handle names

% Save output
function saveDeployment(obj, ev)
f = get(obj,'parent');
h = getappdata(f,'Controls');
doworkspace = ~isempty(regexp(get(obj,'string'),'[Ww]orkspace','once'));
buses = get(h.buses, 'Data');
if(doworkspace)
	% write it to the workspace
	assignin('base','buses',buses);
	msgbox('Successfully saved to workspace');
else
	% save the substructs as .mat files in the current directory
	save([pwd '/buses.mat'], 'buses');
	msgbox(sprintf('Successfully written to file %s', pwd));
end

end

%% Other functions
% size of a degree, in meters
function [latsize, lonsize] = sizeOfADegree(latitude)
latsize = bu.science.astronomy.Earth.lengthOfDegreeLatitude_d( latitude );
lonsize = bu.science.astronomy.Earth.lengthOfDegreeLongitude_d( latitude );
end
