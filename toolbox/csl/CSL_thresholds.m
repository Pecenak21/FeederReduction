function h = CSL_thresholds(imager, time)
%% Create figure elements
h = figure('position', [10,300,1650,550/0.8]);

a(1) = axes('position', [0,0.2,0.333,0.8]);
a(2) = axes('position', [0.3335,0.2,0.333,0.8]);
a(3) = axes('position', [0.667,0.2,0.333,0.8]);
set(a,'ytick',[],'xtick',[]);

c(1) = uicontrol('style', 'slider', 'position', [0, 100, 1650, 20] );
c(2) = uicontrol('style', 'slider', 'position', [0, 60, 1650, 20] );
set(c, 'max', 1.5, 'min', 0, 'sliderstep', [0.0001 0.1]/3, 'value', 0.5, 'callback', @regenplots);
c(3) = uicontrol('String', 'Switch to RBR', 'position', [(1650-100)/2, 20, 100, 20], 'callback', @toggleRBD);


t(1) = uicontrol('style', 'text', 'String', 'Thick Cutoff: 0.5', 'unit', 'pixels', 'Position', [0 122 1650 15]);
t(2) = uicontrol('style', 'text', 'String', 'Thin Cutoff: 0.5', 'unit', 'pixels', 'Position', [0 82 1650 15]);

d(1) = uicontrol('String', '<-1hr', 'position', [(1650-100)/2-300, 20, 80, 20], 'callback', @changeTime);
d(2) = uicontrol('String', '<-10m', 'position', [(1650-100)/2-200, 20, 80, 20], 'callback', @changeTime);
d(3) = uicontrol('String', '<- 1m', 'position', [(1650-100)/2-100, 20, 80, 20], 'callback', @changeTime);
d(4) = uicontrol('String', '1m ->', 'position', [(1650-100)/2+120, 20, 80, 20], 'callback', @changeTime);
d(5) = uicontrol('String', '10m->', 'position', [(1650-100)/2+220, 20, 80, 20], 'callback', @changeTime);
d(6) = uicontrol('String', '1hr->', 'position', [(1650-100)/2+320, 20, 80, 20], 'callback', @changeTime);


%% Set data to the figure
setappdata(h, 'time', time);
setappdata(h, 'imager', imager);
setappdata(h, 'controls', [a c d]);

[img, rbrO, rbrD] = loadnewimage(h);

setappdata(h, 'text', t);

% draw the plots the first time
axes(a(1)); image(img);
axes(a(2)); imagesc(rbrD,[0 1]);
set(a,'DataAspectRatio',size(img));
colorbar
regenplots(c(1));
set(a, 'visible', 'off', 'nextplot', 'replacechildren');
set(a(1), 'nextplot', 'add');

setappdata(h, 'link', linkprop(a,{'XLim','YLim'}));

end

function regenplots(~, ~)

h = gcbf();
if( isempty(h) ) % this just happens when we call it the first time to setup
	h = gcf;
end

c = getappdata(h,'controls');
cutoffs = get(c(4:5),'value');
t = getappdata(h,'text');
set(t(1), 'String', sprintf('Thick Cutoff: %4f',cutoffs{1}));
set(t(2), 'String', sprintf('Thick Cutoff: %4f',cutoffs{2}));

rbrD = getappdata(h,'rbrD');

clouds = (rbrD > cutoffs{2} & rbrD < cutoffs{1}) + 2 * (rbrD > cutoffs{1}); % thick plus thin

% plot the cloudmap it
naturalskymap = [ 0 .4 .8; 1 1 1; .7 .7 .7 ];
z = savezoom(c(3));
axes(c(3));
subimage(uint8(clouds), naturalskymap);
setzoom(c(3),z);
% 

end

function toggleRBD(obj, ~)

h = gcbf();
if( isempty(h) ) % this just happens when we call it the first time to setup
	h = gcf;
end

c = getappdata(h,'controls');
axes(c(2));
z = savezoom(c(2));
s = get(obj,'string');
if(regexp(s,'RBR')) % switch to rbr
	imagesc(getappdata(h,'rbrO'), [0 1.5]);
	set(obj,'string',regexprep(s,'RBR','RBD'));
else % switch to rbd
	imagesc(getappdata(h,'rbrD'), [0 1]);
	set(obj,'string',regexprep(s,'RBD','RBR'));
end
setzoom(c(2),z);

end

function changeTime(obj, ~)

h = gcbf();
if( isempty(h) ) % this just happens when we call it the first time to setup
	h = gcf;
end

% get some data from the figure
s = get(obj, 'string');
time = getappdata(h, 'time');
a = getappdata(h, 'controls');

% figure out how much time to change by
dt = str2double(regexp(s,'\d+','once','match'));
if(regexp(s,'hr'))
	dt = dt*60;
end % otherwise it's in minutes
dt = dt/24/60;
if(any(s=='<'))
	dt = -dt;
end
time = time + dt;
setappdata(h,'time',time);

% load the new time, and plot that image and the rbr
[img rbrO rbrD] = loadnewimage(h);

axes(a(1)); image(img);
axes(a(2));
s = get(a(6),'string');
if(regexp(s,'RBR')) % is rbd (button says 'switch to rbr')
	imagesc(rbrD,[0 1]);
else
	imagesc(rbrO, [0 1.5]);
end

% and regernate the cloudmap
regenplots();

end

function z = savezoom(a)
z = get(a,{'xlim','ylim'});
end
function setzoom(a,z)
set(a,'xlim',z{1},'ylim',z{2});
end

function [img rbrO rbrD] = loadnewimage(h)
imager = getappdata(h,'imager');
time = getappdata(h,'time');

img = imager.readImage(time);
rbrO = rbr(img, imager.darkoffset);
rbrD = imager.cslForTime(time);
rbrD = rbrO - rbrD(:,:,4);
img = previewHDR(img);

setappdata(h, 'image', img);
setappdata(h, 'rbrO', rbrO);
setappdata(h, 'rbrD', rbrD);

end
