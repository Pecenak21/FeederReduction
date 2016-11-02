 function [lineHandles] =  plotSite(feeder)
%% case just a feeder
[pv,bus1,bus2,substation,LinesToPlot,f] = get_info_lines(feeder);
 %% Plot lines
lineHandles = plot([bus1.coord(:,1),bus2.coord(:,1)]',[bus1.coord(:,2),bus2.coord(:,2)]','Color',[0 0 0],'LineWidth',0.1);
%     set(gcf,'Position', [320 -150 1272 624]);
hold on;
set(lineHandles,{'DisplayName'},strcat('line.',{LinesToPlot.name})','LineWidth',1);
set(gcf,'Position',[360 -240 1080 751]);

%% Plot Substation Location
plot(substation(1),substation(2),'v','MarkerSize',20,'MarkerFaceColor','r','MarkerEdgeColor','r');
plot(substation(1),substation(2 ),'^','MarkerSize',20,'MarkerFaceColor','r','MarkerEdgeColor','r');


%% google maps
% plot_google_map('Refresh', 1);hold on;

