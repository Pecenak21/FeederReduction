function[] = PlotOpenDSS_3D_peak(x,y,z,...
    Rotation1,Rotation2,Edge,Label_TimeAxis,Label_DistanceAxis,...
    Label_DataAxis,DataMax,ColorBar,X_Label,Y_Label,Z_Label,Label_FontSize,...
    z_saturation)

%PlotHarmonicData_3D(DataPath,Harmonic_Start,Harmonic_End,Rotation1,Rotation2,Edge,SaveFormat)
%
% PURPOSE : Plots 3D graphs of data from harmonic analysis
%
%
% INPUT :
%
% DataPath : text string with the data path and file name containing the
% data, the file has to be in ASCII format (.txt, .csv,...), the first
% column should contain the frequency and the other columns should contain
% the impedance data
%
% Harmonic_Start : lower limit of the harmonic axis
%
% Harmonic_End : Upper limit of the harmonic axis
%
% Rotation1 : Viewing angle 1, 50 produces a good result
%
% Rotation2 : Viewing angle 2, 30 produces a good result
%
% Edge : Edge = 0 -> no contour line (preferred for a large number of cases)
%        Edge = 1 -> contour line (preferred for a small number of cases)
%
% SaveFormat : 'Print' or 'SaveAs', 
% 'SaveAs' (default) creates .tif and .fig files at the location of the source data
% 'Print' creates .tif and .fig at the location specified in the code
% 'Print' gives you more freedom regarding figure resolution but has
% restrictions regarding data path, e.g., the data path should not contain
% spaces


if nargin == 9
    DataMax='auto';
end

% NumberOfCases = length(Data(1,:));


% CaseNumber=0:(NumberOfCases-1);

if Edge==0
    surf(x,y,z,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
else
    surf(x,y,z,'FaceColor','interp','FaceLighting','phong');
end

% daspect([5 5 1])
axis tight
% view(-50,30)
view(Rotation1,Rotation2)
%camlight left
brighten(jet,.9);
% colorbar
% colormap(jet(128))

if strcmp(Label_TimeAxis,'auto')==0
    set(gca,'ytick',Label_TimeAxis)
end
if strcmp(Label_DistanceAxis,'auto')==0
    set(gca,'xtick',Label_DistanceAxis)
end
if strcmp(Label_DataAxis,'auto')==0
    set(gca,'ztick',Label_DataAxis)
end

if strcmp(DataMax,'auto')==0
    Data_max=DataMax;
%     zlim([0 Data_max]);

    zlim([min(Label_DataAxis) max(Label_DataAxis)]);
    set(gca,'ztick',Label_DataAxis);
%     set(gca,'ZDir','reverse')
%     ColorVector_sample_max=ceil(length(jet(512))*max(max(z))/Data_max);
%     ColorVector_sample_max=ceil(length(jet(512))*max(max(z))/z_saturation);
    ColorVector_sample_max=ceil(length(jet(512))*Data_max/z_saturation(2));
%     ColorVector_sample_max=ceil(length(jet(512))*0.5);
%     ColorVector=fliplr(jet(512)')';
    ColorVector=jet(512);
    ColorVector=ColorVector(1:(ColorVector_sample_max),:);
    colormap(ColorVector);
    freezeColors
else
    colormap(jet(512));
end

if ColorBar==1
%     caxis([min(Label_DataAxis) max(Label_DataAxis)]);
    caxis(z_saturation);
    h=colorbar;
    ax=get(h);
%     keyboard
%     ax(1)=ax(1)+0.5;
%     ax=ax-0.45;
    set(h,'fontsize',9);
%     set(h,'YLim', [0.9 1.2]); 
%     set(h,'Location','East');

    cbfreeze
%     keyboard
end

% caxis([min(Label_DataAxis) max(Label_DataAxis)]);
caxis(z_saturation);
% cbfreeze

set(gca,'FontSize',9)

hx=xlabel(X_Label,'FontSize',Label_FontSize);
set(hx,'Rotation',-17);
hy=ylabel(Y_Label,'FontSize',Label_FontSize);
set(hy,'Rotation',11);
zlabel(Z_Label,'FontSize',Label_FontSize);
% keyboard


% saveas(gcf,[OutputFileName '.tif'])
% saveas(gcf,[OutputFileName '.fig'])