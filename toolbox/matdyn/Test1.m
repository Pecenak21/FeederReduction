function Test1

% Test1
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

%% Modified Euler
% Set options
mdopt(1)=1;     % Modified Euler 

% Run dynamic simulation
fprintf('> Modified Euler...')
[Angles1,Speeds,Eq_tr,Ed_tr,Efd,PM,Voltages,Stepsize1,Errest,Time1]=rundyn('casestagg','casestaggdyn','staggevent',mdopt);

%% Runge-Kutta
% Set options
mdopt(1)=2;     % Runge-Kutta 

% Run dynamic simulation
fprintf('Done.\n> Runge-Kutta...')
[Angles2,Speeds,Eq_tr,Ed_tr,Efd,PM,Voltages,Stepsize2,Errest,Time2]=rundyn('casestagg','casestaggdyn','staggevent',mdopt);

%% Fehlberg
% Set options
mdopt(1)=3;     % Runge-Kutta Fehlberg
mdopt(2)=1e-4;  % tol = 1e-4
mdopt(3)=1e-4;  % minimum step size = 1e-4

% Run dynamic simulation
fprintf('Done.\n> Runge-Kutta Fehlberg...')
[Angles3,Speeds,Eq_tr,Ed_tr,Efd,PM,Voltages,Stepsize3,Errest,Time3]=rundyn('casestagg','casestaggdyn','staggevent',mdopt);
% 
%% Higham-Hall
% Set options
mdopt(1)=4;     % Runge-Kutta Higham-Hall
mdopt(2)=1e-4;  % tol = 1e-4
mdopt(3)=1e-4;  % minimum step size = 1e-4

% Run dynamic simulation
fprintf('Done.\n> Runge-Kutta Higham-Hall...')
[Angles4,Speeds,Eq_tr,Ed_tr,Efd,PM,Voltages,Stepsize4,Errest,Time4]=rundyn('casestagg','casestaggdyn','staggevent',mdopt);
fprintf('Done.\n')

%% Plots
% Plot angles
close all
figure
hold on
xlabel('Time [s]')
ylabel('Generator angles [deg]')
p1 = plot(Time1,Angles1(:,1:2),'-.b');
p2 = plot(Time2,Angles2(:,1:2),':r');
p3 = plot(Time3,Angles3(:,1:2),'--g');
p4 = plot(Time4,Angles4(:,1:2),'m');
Group1 = hggroup;
Group2 = hggroup;
Group3 = hggroup;
Group4 = hggroup;
set(p1,'Parent',Group1)
set(p2,'Parent',Group2)
set(p3,'Parent',Group3)
set(p4,'Parent',Group4)
set(get(get(Group1,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(Group2,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(Group3,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(Group4,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
legend('Modified Euler','Runge-Kutta','Fehlberg','Higham-Hall');
axis([0 Time1(end) -1 1])
axis 'auto y'

figure
hold on
p1 = plot(Time1,Stepsize1,':b');
p2 = plot(Time3,Stepsize3,'--g');
p3 = plot(Time4,Stepsize4,'m');
Group1 = hggroup;
Group2 = hggroup;
Group3 = hggroup;
set(p1,'Parent',Group1)
set(p2,'Parent',Group2)
set(p3,'Parent',Group3)
set(get(get(Group1,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(Group2,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
set(get(get(Group3,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');  
xlabel('Time [s]')
ylabel('Step size')
legend('Modified Euler and Runge-Kutta','Fehlberg','Higham-Hall');
axis([0 Time1(end) -1 1])
axis 'auto y'

return;