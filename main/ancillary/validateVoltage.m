function [fh,o,d] = validateVoltage(c,feederName)
% supported feeder name: fallbrook (520), ramona (971), pointloma (480),
% alpine (355), valleycenter (909)
%

%% plot volt profile for the modelled circuit
[o,d,fh] = plotVoltageProfile(c);

feederNorminalVolt = 12.47; % kv

%% load the voltage data from SDG&E
if ischar(feederName), feederName = lower(feederName); end
switch feederName
    case {'pointloma','480','cabrillo','b',480}
        NumCc='55_0480';
    case {'fallbrook','520','avocado','a',520}
        NumCc='55_0520';
    case {'alpine','355','d',355}
        NumCc='55_0355';
    case {'valleycenter','909','c',909}
        NumCc='54_0909';
    case {'ramona','971','creelman','e',971}
        NumCc='55_0971';
    otherwise
        error('Not supported feeder. Please check or write new code to support this feeder.');
end

x = excel2obj(sprintf('dssconversion/custdata/%s_-_Balanced_Results.csv', NumCc));
fn=fieldnames(x);x=x.(fn{1}); 

%% Calculate Voltage
S = sqrt([x.Into_kW].*[x.Into_kW]+[x.Into_kvar].*[x.Into_kvar]); % kVA
I = [x.Into_Amps]; % Amps
Dist = [x.Section_Dist]/3.28084; % km
V_LL = S./I/sqrt(3)/feederNorminalVolt; % Line to neutral voltage normalized by feeder norminal voltage
aa = horzcat(Dist',V_LL');
aa = aa(~isnan(aa(:,2)),:);

%% pick only voltage in the .85-1.15 range for now
% 12kV range 3 phases
a12 = aa((aa(:,2)>0.85) & (aa(:,2)<1.15),:);
% 12kV range phase neutral
% a67 = aa((aa(:,2)>12/sqrt(3)*0.9) & (aa(:,2)<12/sqrt(3)*1.1),:);
% % 4.16kV range
% a416 = aa((aa(:,2)>4.16*0.9) & (aa(:,2)<4.16*1.1),:);
% % 4.16kV range phase neutral
% a24 = aa((aa(:,2)>4.16/sqrt(3)*0.9) & (aa(:,2)<4.16/sqrt(3)*1.1),:);

%% Plot
[~, order] = sort(a12(:,1));
figure(fh); plot(a12(order,1),a12(order,2),'.k');
legend({'Phase 1','Phase 2','Phase 3','Utility data'});

end