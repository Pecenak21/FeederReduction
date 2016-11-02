function [ figHandle ] = resetFigure( figHandle )
% resetFigure for publishing purpose
set(figHandle,'units', 'pixels');
set(figHandle,'position', [420   270   440   309]);
set(figHandle,'color','w');
grid on; box on;
end

