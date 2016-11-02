function [ x ] = plotFeederProfileSynergee( FeederNum )
%VOLTAGEDROP Summary of this function goes here
%   Detailed explanation goes here
%   INPUT:
%       - "FeederNum": Number of the feeder (355, 480, 520...)
%       - "outpath": Path to the file where to save the chart

%% load the right Balanced Results
switch FeederNum
    case 480
        NumCc='55_0480';
    case 520
        NumCc='55_0520';
    case 355
        NumCc='55_0355';
    case 909
        NumCc='54_0909';
    case 971
        NumCc='55_0971';
end

x = excel2obj(sprintf('dssconversion/custdata/%s_-_Balanced_Results.csv', NumCc));
fn=fieldnames(x);x=x.(fn{1}); 

%% Voltage profile
figure, plot([x.Section_Dist]*.3048,[x.Volts_Out]/120,'.');
grid on; box on; set(gcf,'color','w'); %ylim([.9 1.1]);
set(gca,'fontsize',12);
title('SYNERGEE Voltage profile'); xlabel('Distance [km]'); ylabel('Voltage [p.u.]'); 
set(gca,'xtick',[0:2:20]); set(gca,'ylim',[.8 1.2]);
%% current profile
figure, plot([x.Section_Dist]*.3048,[x.Into_Amps],'.');
grid on; box on; set(gcf,'color','w'); %ylim([.9 1.1]);
set(gca,'fontsize',12);
title('SYNERGEE Current profile'); xlabel('Distance [km]'); ylabel('Current [Amps]'); 
set(gca,'xtick',[0:2:20]);

%% Power factor
figure, plot([x.Section_Dist]*.3048,[x.Into_pf]/100,'.'); hold on;
id = [x.Into_pf] <0;
plot([x(id).Section_Dist]*.3048,abs([x(id).Into_pf]/100),'.r');
grid on; box on; set(gcf,'color','w'); %ylim([.9 1.1]);
set(gca,'fontsize',12);
title('SYNERGEE Power factor'); xlabel('Distance [km]'); ylabel('Power Factor [-]'); 
set(gca,'xtick',[0:2:20]);

%% total power
figure, plot(	[x.Section_Dist]*.3048,sqrt([x.Into_kvar].^2 + [x.Into_kW].^2),'.',...
				[x.Section_Dist]*.3048,[x.Into_kW],'.',...
				[x.Section_Dist]*.3048,[x.Into_kvar],'.');
legend({'Total Apparent Power [kVA]','Active Power [kW]','Reactive Power [kVar]'});
grid on; box on; set(gcf,'color','w'); %ylim([.9 1.1]);
set(gca,'fontsize',12);
title('SYNERGEE Power'); xlabel('Distance [km]'); ylabel('Power [kX]'); 
set(gca,'xtick',[0:2:20]);

end
