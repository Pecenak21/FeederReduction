function Test2

% Test2
% 
% test script

% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%%
mdopt=Mdoption;
mdopt(5)=0;     % no progress info
mdopt(6)=0;     % no plots

%% Fehlberg 1
% Set options
mdopt(1)=3;     % Runge-Kutta Fehlberg
mdopt(2)=1e-4;  % tol = 1e-4
mdopt(3)=1e-4;  % minimum step size = 1e-4

% Run dynamic simulation
fprintf('> Runge-Kutta Fehlberg...')
[Angles1,Speeds,Eq_tr,Ed_tr,Efd,PM,Voltages,Stepsize1,Errest,Time1]=rundyn('casestagg','casestaggdyn','staggevent',mdopt);
% 
%% Fehlberg 2
% Set options
mdopt(1)=3;     % Runge-Kutta Fehlberg
mdopt(2)=1e-3;  % tol = 1e-3
mdopt(3)=1e-4;  % minimum step size = 1e-4

% Run dynamic simulation
fprintf('Done.\n> Runge-Kutta Fehlberg...')
[Angles2,Speeds,Eq_tr,Ed_tr,Efd,PM,Voltages,Stepsize2,Errest,Time2]=rundyn('casestagg','casestaggdyn','staggevent',mdopt);
fprintf('Done.\n')

%% Plots
% Plot angles
% close all
figure
hold on
xlabel('Time [s]')
ylabel('Generator angles [deg]')
p1 = plot(Time1,Angles1,'b');
p2 = plot(Time2,Angles2,'r--');
Group1 = hggroup;
Group2 = hggroup;
set(p1,'Parent',Group1)
set(p2,'Parent',Group2)
set(get(get(Group1,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(Group2,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
legend('RK Fehlberg 1e-4','RK Fehlberg 1e-3');
axis([0 Time1(end) -1 1])
axis 'auto y'

figure
hold on
xlabel('Time [s]')
ylabel('Step size')
p1 = plot(Time1,Stepsize1,'b');
p2 = plot(Time2,Stepsize2,'r--');
Group1 = hggroup;
Group2 = hggroup;
set(p1,'Parent',Group1)
set(p2,'Parent',Group2)
set(get(get(Group1,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(Group2,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
legend('RK Fehlberg 1e-4','RK Fehlberg 1e-3');
axis([0 Time1(end) -1 1])
axis 'auto y'

return;