clear all
close all
clc

%exctracts data from the figures obtained from the daily simulation runs (DailySimulations.m)
%calculates the average voltage margin between case with pv and case
%without pv over the period of time with daylight (7am - 4pm).


strDirectory='C:\Work\Projects\2012\1787-UCSD_PV\Simulation\System\SDGE';
cd(strDirectory); 
curdir = pwd;
for i=1:4
    switch i
        case 1
            path = [curdir '/tmp_scnr1'];
            for j=0:6
                disp(['Scenario ' num2str(i), ' Case ' num2str(j)]);
                %Extract voltages
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Maximum and Minimum Bus Voltages.fig']);
                D = get(gca,'Children');
                Child = get(D);
                %extract data from figure
                MaxPV = Child(4).YData;MaxNoPV = Child(3).YData;MinPV = Child(2).YData;MinNoPV = Child(1).YData;
                min=7;max=16; %limit the data between 7am and 4pm (sunlight hours)
                if length(MaxPV) > 24
                    min=840;max=1920;%limits for 30sec resolution data
                end
                AverageDivMax = sum(MaxPV(min:max)-MaxNoPV(min:max))/length(MaxPV(min:max))*100;
                AverageDivMin = sum(MinPV(min:max)-MinNoPV(min:max))/length(MinPV(min:max))*100;
                disp([' (max(V)): ' num2str(AverageDivMax) ' and (min(V)): ' num2str(AverageDivMin)]);
                close all;
               
                %Extract losses
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Losses.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                h=1;%hourly devider
                if length(MaxPV) > 24
                    h=120;%hourly devider for 30sec resolution data
                end
                MWhNoPV = sum(MWNoPV)/h;MWhPV=sum(MWPV)/h;MVArhNoPV=sum(MVarNoPV)/h;MVArhPV=sum(MVarPV)/h;
                disp(['daily losses: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                MWhNoPV = sum(MWNoPV(min:max))/h;MWhPV=sum(MWPV(min:max))/h;MVArhNoPV=sum(MVarNoPV(min:max))/h;MVArhPV=sum(MVarPV(min:max))/h;
                disp(['daylight losses: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                close all;
               
                %Extract power
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Power.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                MWhNoPV = sum(MWNoPV)/h;MWhPV=sum(MWPV)/h;MVArhNoPV=sum(MVarNoPV)/h;MVArhPV=sum(MVarPV)/h;
                disp(['daily consumption: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                MWhNoPV = sum(MWNoPV(min:max))/h;MWhPV=sum(MWPV(min:max))/h;MVArhNoPV=sum(MVarNoPV(min:max))/h;MVArhPV=sum(MVarPV(min:max))/h;
                disp(['daylight consumption: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                close all;
            end
        case 2
            disp('________________________________________________________________');
            path = [curdir '/tmp_scnr2'];
            for j=0:6
                disp(['Scenario ' num2str(i), ' Case ' num2str(j)]);
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Maximum and Minimum Bus Voltages.fig']);
                D = get(gca,'Children');
                Child = get(D);
                %extract data from figure
                MaxPV = Child(4).YData;MaxNoPV = Child(3).YData;MinPV = Child(2).YData;MinNoPV = Child(1).YData;
                min=7;max=16; %limit the data between 7am and 4pm (sunlight hours)
                if length(MaxPV) > 24
                    min=840;max=1920;%limits for 30sec resolution data
                end
                AverageDivMax = sum(MaxPV(min:max)-MaxNoPV(min:max))/length(MaxPV(min:max))*100;
                AverageDivMin = sum(MinPV(min:max)-MinNoPV(min:max))/length(MinPV(min:max))*100;
                disp([' (max(V)): ' num2str(AverageDivMax) ' and (min(V)): ' num2str(AverageDivMin)]);
                close all;
               
                %Extract losses
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Losses.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                h=1;%hourly devider
                if length(MaxPV) > 24
                    h=120;%hourly devider for 30sec resolution data
                end
                MWhNoPV = sum(MWNoPV)/h;MWhPV=sum(MWPV)/h;MVArhNoPV=sum(MVarNoPV)/h;MVArhPV=sum(MVarPV)/h;
                disp(['daily losses: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                MWhNoPV = sum(MWNoPV(min:max))/h;MWhPV=sum(MWPV(min:max))/h;MVArhNoPV=sum(MVarNoPV(min:max))/h;MVArhPV=sum(MVarPV(min:max))/h;
                disp(['daylight losses: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                close all;
               
                %Extract power
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Power.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                MWhNoPV = sum(MWNoPV)/h;MWhPV=sum(MWPV)/h;MVArhNoPV=sum(MVarNoPV)/h;MVArhPV=sum(MVarPV)/h;
                disp(['daily consumption: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                MWhNoPV = sum(MWNoPV(min:max))/h;MWhPV=sum(MWPV(min:max))/h;MVArhNoPV=sum(MVarNoPV(min:max))/h;MVArhPV=sum(MVarPV(min:max))/h;
                disp(['daylight consumption: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                close all;
            end
        case 3
            disp('________________________________________________________________');
            path = [curdir '/tmp_scnr3'];
            for j=0:6
                disp(['Scenario ' num2str(i), ' Case ' num2str(j)]);
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Maximum and Minimum Bus Voltages.fig']);
                D = get(gca,'Children');
                Child = get(D);
                %extract data from figure
                MaxPV = Child(4).YData;MaxNoPV = Child(3).YData;MinPV = Child(2).YData;MinNoPV = Child(1).YData;
                min=7;max=16; %limit the data between 7am and 4pm (sunlight hours)
                if length(MaxPV) > 24
                    min=840;max=1920;%limits for 30sec resolution data
                end
                AverageDivMax = sum(MaxPV(min:max)-MaxNoPV(min:max))/length(MaxPV(min:max))*100;
                AverageDivMin = sum(MinPV(min:max)-MinNoPV(min:max))/length(MinPV(min:max))*100;
                disp([' (max(V)): ' num2str(AverageDivMax) ' and (min(V)): ' num2str(AverageDivMin)]);
                close all;
               
                %Extract losses
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Losses.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                h=1;%hourly devider
                if length(MaxPV) > 24
                    h=120;%hourly devider for 30sec resolution data
                end
                MWhNoPV = sum(MWNoPV)/h;MWhPV=sum(MWPV)/h;MVArhNoPV=sum(MVarNoPV)/h;MVArhPV=sum(MVarPV)/h;
                disp(['daily losses: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                MWhNoPV = sum(MWNoPV(min:max))/h;MWhPV=sum(MWPV(min:max))/h;MVArhNoPV=sum(MVarNoPV(min:max))/h;MVArhPV=sum(MVarPV(min:max))/h;
                disp(['daylight losses: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                close all;
               
                %Extract power
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Power.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                MWhNoPV = sum(MWNoPV)/h;MWhPV=sum(MWPV)/h;MVArhNoPV=sum(MVarNoPV)/h;MVArhPV=sum(MVarPV)/h;
                disp(['daily consumption: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                MWhNoPV = sum(MWNoPV(min:max))/h;MWhPV=sum(MWPV(min:max))/h;MVArhNoPV=sum(MVarNoPV(min:max))/h;MVArhPV=sum(MVarPV(min:max))/h;
                disp(['daylight consumption: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                close all;
            end
       case 4
            disp('________________________________________________________________');
            path = [curdir '/tmp_scnr4'];
            for j=0:6
                disp(['Scenario ' num2str(i), ' Case ' num2str(j)]);
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Maximum and Minimum Bus Voltages.fig']);
                D = get(gca,'Children');
                Child = get(D);
                %extract data from figure
                MaxPV = Child(4).YData;MaxNoPV = Child(3).YData;MinPV = Child(2).YData;MinNoPV = Child(1).YData;
                min=7;max=16; %limit the data between 7am and 4pm (sunlight hours)
                if length(MaxPV) > 24
                    min=840;max=1920;%limits for 30sec resolution data
                end
                AverageDivMax = sum(MaxPV(min:max)-MaxNoPV(min:max))/length(MaxPV(min:max))*100;
                AverageDivMin = sum(MinPV(min:max)-MinNoPV(min:max))/length(MinPV(min:max))*100;
                disp([' (max(V)): ' num2str(AverageDivMax) ' and (min(V)): ' num2str(AverageDivMin)]);
                close all;
               
                %Extract losses
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Losses.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                h=1;%hourly devider
                if length(MaxPV) > 24
                    h=120;%hourly devider for 30sec resolution data
                end
                MWhNoPV = sum(MWNoPV)/h;MWhPV=sum(MWPV)/h;MVArhNoPV=sum(MVarNoPV)/h;MVArhPV=sum(MVarPV)/h;
                disp(['daily losses: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                MWhNoPV = sum(MWNoPV(min:max))/h;MWhPV=sum(MWPV(min:max))/h;MVArhNoPV=sum(MVarNoPV(min:max))/h;MVArhPV=sum(MVarPV(min:max))/h;
                disp(['daylight losses: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                close all;
               
                %Extract power
                open([path '\f520_case' num2str(j) '\results/Case ' num2str(j) '_Total Power.fig']);
                D = get(gca,'Children');
                Child = get(D);
                MVarNoPV=Child(1).YData;MVarPV=Child(2).YData;MWNoPV=Child(3).YData;MWPV=Child(4).YData;
                MWhNoPV = sum(MWNoPV)/h;MWhPV=sum(MWPV)/h;MVArhNoPV=sum(MVarNoPV)/h;MVArhPV=sum(MVarPV)/h;
                disp(['daily consumption: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                MWhNoPV = sum(MWNoPV(min:max))/h;MWhPV=sum(MWPV(min:max))/h;MVArhNoPV=sum(MVarNoPV(min:max))/h;MVArhPV=sum(MVarPV(min:max))/h;
                disp(['daylight consumption: (MWhNoPV): ' num2str(MWhNoPV) ' (MVArhNoPV): ' num2str(MVArhNoPV) ' (MWhPV): ' num2str(MWhPV) ' (MVArhPV): ' num2str(MVArhPV)]);
                close all;
            end
    end
end