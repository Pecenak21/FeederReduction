function drawSections(sectionIds,colorarg,figureid)
% Draws lines for the listed sections in a new figure window
%	You must call this function once with different args to cache the
%	needed data structures:
%	 drawSections('generate',d);
%	Then you can call it repeatedly to graph sections
%
%	 drawSections(sectionIds,colorarg,figureid);
%	
%	figureid (optional) specifies a figure to draw to (default, new figure)
%	colorarg (optional) is passed to plot() as a linestyle argument, so can be any of
%	the usual color strings
%	sectionIDs should be a cell array of strings
%
%	if you specify a colormap name instead of a plot color (i.e. more than
%	2 characters), segment directions are colorized instead.  Complicated
%	color schemes like Jet are NOT recommended (I like winter)
persistent nodeXY;
persistent sectionH;

%% generate the persistent caches if requested
if(strcmp(sectionIds,'generate'))
	d = colorarg;
	%% build a hash lookup table for node coords
	for i=1:length(d.Node)
		nodeXY.(['n' d.Node(i).NodeId]) = [d.Node(i).X d.Node(i).Y];
	end

	%% build a hash lookup table for sections
	for i=1:length(d.InstSection)
		id = ['s' d.InstSection(i).SectionId];
		id(id==' ') = [];
		sectionH.(id) = struct('FromNodeId',d.InstSection(i).FromNodeId,'ToNodeId',d.InstSection(i).ToNodeId);
	end
	return;
end

%% process inputs
if(nargin < 3 || isempty(figureid))
	figure();
else
	figure(figureid);
end
if(nargin < 2 || isempty(colorarg))
	colorarg = 'b';
end
sectionIds = unique(sectionIds);

%% do the drawing
oldop = get(gca,'next');
set(gca,'next','add');
for i=1:length(sectionIds)
	id = ['s' sectionIds{i}];
	id(id==' ') = [];
	section = sectionH.(id);
	coords = vertcat(nodeXY.(['n' section.FromNodeId]),nodeXY.(['n' section.ToNodeId]));
	if(length(colorarg) > 2)
		coords = vertcat(coords,[NaN NaN]);
		id = patch(coords(:,1),coords(:,2),[0 1 1]);
		set(id,'cdata',[1 0 0],'edgecolor','interp');
	else
		plot(coords(:,1),coords(:,2),colorarg);
	end
end
if(length(colorarg) > 2)
	colormap(colorarg);
end
% restore old drawing mode
set(gca,'next',oldop);

end
