function drawGHIFigures(cdata1, fig_title, figureFileName)
%CREATEFIGURE(CDATA1)
%  CDATA1:  image cdata

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YDir','reverse','Layer','top',...
    'DataAspectRatio',[1 1 1]);
xlim([0.5 1002]);
ylim([0.5 1002]);
box('on');
hold('all');

% Create image
thisfig = image(cdata1,'Parent',axes1,'CDataMapping','scaled');
% Create title
%title({'GHI 2009 10 04 hour ', hours, 'minute ', minutes, 'seconds ', seconds, ' '});
title({fig_title});

% write to file
saveas(thisfig, figureFileName,'png');





