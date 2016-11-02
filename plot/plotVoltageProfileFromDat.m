function plotVoltageProfileFromDat( FeederNum, outpath )
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

%% Calculate Voltage
S = sqrt([x.Into_kW].*[x.Into_kW]+[x.Into_kvar].*[x.Into_kvar]); % kVA
I = [x.Into_Amps]; % Amps
Dist = [x.Section_Dist]; % kft
V_LL = S./I/sqrt(3);
aa = horzcat(Dist',V_LL');
aa = aa(~isnan(aa(:,2)),:);

%% sort Voltage ranges
% 12kV range 3 phases
a12 = aa((aa(:,2)>12*0.9) & (aa(:,2)<12*1.1),:);
% 12kV range phase neutral
a67 = aa((aa(:,2)>12/sqrt(3)*0.9) & (aa(:,2)<12/sqrt(3)*1.1),:);
% 4.16kV range
a416 = aa((aa(:,2)>4.16*0.9) & (aa(:,2)<4.16*1.1),:);
% 4.16kV range phase neutral
a24 = aa((aa(:,2)>4.16/sqrt(3)*0.9) & (aa(:,2)<4.16/sqrt(3)*1.1),:);

%% 

%% Plot
[order order] = sort(a12(:,1));
f=figure; plot(a12(order,1)*0.3048,a12(order,2)/12,'.','MarkerSize',10);
title('Voltage function of distance from substation')
xlabel('Distance [km]');
ylabel('Voltage [p.u.]');
grid on;
set(gcf,'color','w');
if ~exist('outpath','var') || isempty(outpath)
	outpath = pwd;
end
saveas(f,[[outpath '/VoltageDist_' NumCc], '.png'])
saveas(f,[[outpath '/VoltageDist_' NumCc], '.fig'], 'fig')
end
