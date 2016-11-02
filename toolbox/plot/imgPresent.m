function [h, a]  = imgPresent(img)
% IMGPRESENT(img) draws the planes of an image in three axes in the same window
%
% This can be used for viewing the RGB components separately, or other
% representations such as HSV or LAB color.
%
% The three axes have their x and y limits linked so that they pan/zoom
% together.  Figure and axes handles can be returned as outputs (first and
% second respectively) for further tweaking.

%% draw it
% create a figure and give it a better position
h_ = figure;
p = get(h_,'position');
dh = 380-p(4);
p(4) = 380;
p(2) = p(2)-dh;
p(3) = 1200;
set(h_,'position',p);

% create three new axes that each fill 1/3 of the window (except for adding their colorbars)
a(1) = axes('position',[0 0 1/3 1]);
imagesc(img(:,:,1)); colorbar;
a(2) = axes('position',[1/3 0 1/3 1]);
imagesc(img(:,:,2)); colorbar;
a(3) = axes('position',[2/3 0 1/3 1]);
imagesc(img(:,:,3)); colorbar;
set(a, 'XTick', [], 'YTick', [], 'DataAspectRatio', [1 1 1]);

% create a new linkprop object and save a copy of it with the figure so that it will continue to link the axes properties until the figure is closed
setappdata(h_,'proplink',linkprop(a,{'xlim','ylim'}));

%% return a figure handle if requested
if(nargout > 0)
	h = h_;
end

end