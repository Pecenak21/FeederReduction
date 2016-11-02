function bg = getBiograph(c)
global indent; if isempty(indent), indent = ''; end
global conf; global fName;
fp = sprintf('%s/%s%s.mat',conf.outputDir,'Biograph_',fName);
if exist(fp,'file')
    fprintf(['%sSaved file for neighbor nodes exists: ' fp '. Load to use.\n'],indent);
    bg = load(fp); bg = bg.dat;
    return;
end

bus = unique([cleanBus({c.line.bus1}) cleanBus({c.line.bus2})]);
% edges are lines
% build connected map based on line infor
conMat = zeros(length(bus),length(bus));
for i = 1:length(c.line)
	[~, id] = ismember(cleanBus({c.line(i).bus1, c.line(i).bus2}),bus);
 	conMat(id(1),id(2)) = c.line(i).Length;
 	conMat(id(2),id(1)) = c.line(i).Length;
end
for i = 1:length(c.transformer)
	if ismember('sourcebus',lower(cleanBus(c.transformer(i).bus)))
		continue;
	end
	[~, id] = ismember(cleanBus(c.transformer(i).bus),bus);
 	conMat(id(1),id(2)) = 1;
 	conMat(id(2),id(1)) = 1;
end
bg = biograph(conMat,bus,'ShowWeights','on');
% bg.ShowArrows = 'off';
% bg.EdgeType = 'straight';
%view(bg)
% find neighbors of Generation/PV nodes (that also have storage systems)
saveFile(fp,bg);
fprintf(['%sSaved neighbor nodes'' file: ' fp '\n'],indent);
end