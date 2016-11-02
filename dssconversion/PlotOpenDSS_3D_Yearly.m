function[] = PlotOpenDSS_3D(Values_DistanceAxis,... 
Values_DataAxis_NoPV,...
Values_DataAxis_WithPV,... 
DirectoryOutput,...
TimeMin,...
TimeMax,...
Edge,...
X_Label_2D,...
Y_Label_2D,...
Z_Label_2D,...
X_Label_3D,...
Y_Label_3D,...
Z_Label_3D,...
Label_FontSize_2D,...
Label_FontSize_3D,...
z_saturation,...
z_range,...
SimName,...
nth,...
xrange,...
DTick)


Title = {'No PV', 'With PV'};

disp(['- Plotting ' SimName]);

Directory_Base_Output = DirectoryOutput;
PlotName = SimName;
DataPath_Output = [Directory_Base_Output '/' PlotName];
% y_time=(TimeMin:nth:TimeMax);
z_NoPV = Values_DataAxis_NoPV;
z_WithPV = Values_DataAxis_WithPV;
    
z_data_total = struct('z_NoPV_data',z_NoPV,'z_WithPV_data',z_WithPV);

OutputFileName = DataPath_Output;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjusting time axis labeling
% Label_TimeAxis=(TimeMin:TimeMax/8:TimeMax); % labeling of the time axis
if nth>1
%     y_time = y_time./120;
%     Label_TimeAxis = round(Label_TimeAxis./120);
    y_time=xrange(1:nth:end);
    time_step = round(length(xrange)/5);
    Label_TimeAxis = xrange(1:time_step:end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjusting distance axis labeling
RoundingFactor_x=10;
x_max=max(Values_DistanceAxis);
x_distance=Values_DistanceAxis;
step_x=ceil(x_max/RoundingFactor_x)*RoundingFactor_x/10;
Label_DistanceAxis_3d=0:step_x:x_max;% labeling of the 3d axis
Label_DistanceAxis_2d=0:step_x:x_max;% labeling of the 2d axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjusting data axis labeling
% z_range = [0.85:0.04:1.05];
Label_DataAxis=z_range;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z_1 = [z_data_total.z_NoPV_data; z_data_total.z_WithPV_data;];
b=0;c=0;
f=figure('visible','on');
       
for a=1:2
    z_data = z_1(a+b:TimeMax/nth*a,:);
    z_max=z_saturation(2);    
% 3D plot
    h = subplot(4,2,a);
    ax=get(h,'Position'); 
    ax(1)=ax(1)-(0.06+c);ax(2)=ax(2)-0.03;ax(3)=ax(3)+0.08;ax(4)=ax(4)+0.03; % adjusts subplot location within the figure
    set(h,'Position',ax);hold on
    PlotOpenDSS_3D_peak(x_distance,y_time,z_data,50,30,Edge,Label_TimeAxis,Label_DistanceAxis_3d,Label_DataAxis,z_max,1,X_Label_3D,Y_Label_3D,Z_Label_3D,Label_FontSize_3D,z_saturation)
    PlotOpenDSS_3D_peak(x_distance,y_time,z_data,50,30,Edge,Label_TimeAxis,Label_DistanceAxis_3d,Label_DataAxis,z_max,0,X_Label_3D,Y_Label_3D,Z_Label_3D,Label_FontSize_3D,z_saturation)
    grid on
    title(char(Title(a)), 'FontWeight', 'bold', 'FontSize', 12);
    datetick('y',DTick,'keeplimits', 'keepticks')
    hold off

% Voltage vs. Time plot (view from right)
    time_step = round(length(xrange)/6);
    Label_TimeAxis = xrange(1:time_step:end);% add more ticks to the time axis - just looks better
    h = subplot(4,2,a+2); 
    ax=get(h,'Position'); 
    ax(1)=ax(1)-(0.06+c);ax(2)=ax(2)-0.045;ax(3)=ax(3)+0.03;ax(4)=ax(4)+0.01;  
    set(h,'Position',ax);hold on
    PlotOpenDSS_2D_peak(x_distance,y_time,z_data,90,0,Label_TimeAxis,Label_DistanceAxis_2d,Label_DataAxis,z_max,0,X_Label_2D,Y_Label_2D,Z_Label_2D,Label_FontSize_2D,z_saturation)
    PlotOpenDSS_2D_peak(x_distance,y_time,z_data,90,0,Label_TimeAxis,Label_DistanceAxis_2d,Label_DataAxis,z_max,0,X_Label_2D,Y_Label_2D,Z_Label_2D,Label_FontSize_2D,z_saturation)
    grid on
    datetick('y',DTick,'keeplimits', 'keepticks')
    hold off
    
% Time vs. Distance (view from top)
    h = subplot(4,2,a+4); 
    ax=get(h,'Position'); 
    ax(1)=ax(1)-(0.06+c);ax(2)=ax(2)-0.055;ax(3)=ax(3)+0.03;ax(4)=ax(4)+0.01;  
    set(h,'Position',ax);hold on
    PlotOpenDSS_2D_peak(x_distance,y_time,z_data,0,90,Label_TimeAxis,Label_DistanceAxis_2d,Label_DataAxis,z_max,0,X_Label_2D,Y_Label_2D,Z_Label_2D,Label_FontSize_2D,z_saturation)
    PlotOpenDSS_2D_peak(x_distance,y_time,z_data,0,90,Label_TimeAxis,Label_DistanceAxis_2d,Label_DataAxis,z_max,0,X_Label_2D,Y_Label_2D,Z_Label_2D,Label_FontSize_2D,z_saturation)
    grid on
    datetick('y',DTick,'keeplimits', 'keepticks')
    axis tight
    hold off
      
% Voltage vs. Distance (view from front)
    h = subplot(4,2,a+6);
    ax=get(h,'Position'); 
    ax(1)=ax(1)-(0.06+c);ax(2)=ax(2)-0.065;ax(3)=ax(3)+0.03;ax(4)=ax(4)+0.01;  
    set(h,'Position',ax);hold on
    PlotOpenDSS_2D_peak(x_distance,y_time,z_data,0,0.01,Label_TimeAxis,Label_DistanceAxis_2d,Label_DataAxis,z_max,0,X_Label_2D,Y_Label_2D,Z_Label_2D,Label_FontSize_2D,z_saturation)
    PlotOpenDSS_2D_peak(x_distance,y_time,z_data,0,0.01,Label_TimeAxis,Label_DistanceAxis_2d,Label_DataAxis,z_max,0,X_Label_2D,Y_Label_2D,Z_Label_2D,Label_FontSize_2D,z_saturation)
    grid on
    hold off
     
    b=b+TimeMax/nth-1;c=c-0.05;
end

subtitle(PlotName); % plots the main title on the figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjusts the figure to fit the computer screen
screen=get(0,'screensize');
offset1=34;
offset2=111;
set(gcf,'units','pixels','Position',[2 offset1 screen(3) screen(4)-offset2])
set(gcf, 'PaperPositionMode', 'auto')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveas(f,[OutputFileName '-3D.fig'], 'fig');
print('-djpeg', '-noui', '-r0', [OutputFileName '-3D.png']); % saves the figure
close all
end


