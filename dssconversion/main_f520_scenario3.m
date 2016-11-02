close all
clear all
pack
clc

% strDirectory='C:\Work\Projects\2012\1787-UCSD_PV\Simulation\System\SDGE';
% cd(strDirectory); 
curdir = pwd;

temp = 'tmp_scnr3';

if ~exist(temp,'dir')
    mkdir(temp);
end

%% Scenario 3: 10 times more PV than scenario 1 (with 2 big pv sites)
%% Case 0: 4 different loadshapes for all the loads and a single loadshape for all the PVs in the system. Time interval: 1 hour 
%% loading ...
clear; 
temp = 'tmp_scnr3';
load('data/f520.mat'); load('data/load.mat');

%% Create load profiles 
%aggregated loadshapes
l1 = la;
[x idx] = ismember({l1.Name}',strcat('Load.',{c.load.Name}'));

dt = 1/24;
t = dt:dt:1;
offset = 0;%length(c.loadshape);
for i = 1:length(l1)
    % get id of loadshape considering the offset
    id = offset + i;
    c.loadshape(id) = dssloadshape;
    c.loadshape(id).Name = ['loadshape_' c.load(idx(i)).Name];
    c.loadshape(id).interval = 1;
    c.loadshape(id).Npts = length(t);
    % set up multiplier; assume the load profile starts at 12am in the morning
    m = [l1(i).x1_00_AM l1(i).x2_00_AM l1(i).x3_00_AM l1(i).x4_00_AM l1(i).x5_00_AM...
        l1(i).x6_00_AM l1(i).x7_00_AM l1(i).x8_00_AM l1(i).x9_00_AM l1(i).x10_00_AM l1(i).x11_00_AM l1(i).x12_00_N...
        l1(i).x1_00_PM l1(i).x2_00_PM l1(i).x3_00_PM l1(i).x4_00_PM l1(i).x5_00_PM l1(i).x6_00_PM...
        l1(i).x7_00_PM l1(i).x8_00_PM l1(i).x9_00_PM l1(i).x10_00_PM l1(i).x11_00_PM l1(i).x12_00_AM];
    % interpolate to get 30 second interval data
    c.loadshape(id).Mult = m;
    c.load(idx(i)).Daily = c.loadshape(id).Name;
end

%% apply pv profile of the 1MW PV site at hospital to all PV systems
pv = load('Fallbrook_Scenario2/scenario2_pvProfiles.mat','pv');pv = pv.pv;
[x pid] = max([pv.kVA]);

dt = 1/24;
t = dt:dt:1;
% create one more PV loadshape
idx = length(c.load) + 1;
c.loadshape(idx) = dssloadshape;
nm = 'pv';
c.loadshape(idx).Name = nm;
c.loadshape(idx).Npts = length(t);
c.loadshape(idx).interval = 1;
m = zeros(1,length(t));
% find data that we have in PV profile and assign it to m taking into account the time zone for California and the collected data date of Dec 14th 2012 
% special handling needed to avoid rounding error for matlab's time/date
[x, id2] = ismember(round((t+ datenum('2012-12-14 00:00:00'))*24), round((pv(pid).time-8/24)*24) );
id2 = id2(id2>0);
m(x) = pv(pid).gi(id2)/1000;
% quality control: set NaN values to 0 to avoid fricking out OpenDSS
m(isnan(m)) = 0;
c.loadshape(idx).Mult = m;

% match PV systems in the circuit with calculated profiles
% load new pv system
pv = load('data/scenario2_pv.mat'); 
c.pvsystem = pv.pv;
for i = 1:length(c.pvsystem)
    % assign newly created profile to according PV system
    c.pvsystem(i).daily = nm;
end

%% END case 0
p = dsswrite(c,'f520',1,[temp '/f520_case0']);
save([temp '/scenario3_case0.mat'],'c','glc');

%% Case 1: 4 different loadshapes for all the loads and a single loadshape for all the PVs in the system. Time interval: 30s 
%% loading ...
clear;
temp = 'tmp_scnr3';
load('data/f520.mat'); load('data/load.mat');

%% Create load profiles 
%aggregated loadshapes
l1 = la;
[x idx] = ismember({l1.Name}',strcat('Load.',{c.load.Name}'));

dt = 30/3600/24;
t = dt:dt:1;
offset = 0;%length(c.loadshape);
for i = 1:length(l1)
    % get id of loadshape considering the offset
    id = offset + i;
    c.loadshape(id) = dssloadshape;
    c.loadshape(id).Name = ['loadshape_' c.load(idx(i)).Name];
    c.loadshape(id).sInterval = 30;
    c.loadshape(id).Npts = length(t);
    % set up multiplier; assume the load profile starts at 12am in the morning
    m = [l1(i).x12_00_AM l1(i).x1_00_AM l1(i).x2_00_AM l1(i).x3_00_AM l1(i).x4_00_AM l1(i).x5_00_AM...
        l1(i).x6_00_AM l1(i).x7_00_AM l1(i).x8_00_AM l1(i).x9_00_AM l1(i).x10_00_AM l1(i).x11_00_AM l1(i).x12_00_N...
        l1(i).x1_00_PM l1(i).x2_00_PM l1(i).x3_00_PM l1(i).x4_00_PM l1(i).x5_00_PM l1(i).x6_00_PM...
        l1(i).x7_00_PM l1(i).x8_00_PM l1(i).x9_00_PM l1(i).x10_00_PM l1(i).x11_00_PM l1(i).x12_00_AM];
    % interpolate to get 30 second interval data
    c.loadshape(id).Mult = interp1(0:1/24:1,m,t,'cubic');
    c.load(idx(i)).Daily = c.loadshape(id).Name;
end

%% apply pv profile of the 1MW PV site at hospital to all PV systems
pv = load('dssconversion/Fallbrook_Scenario2/scenario2_pvProfiles.mat','pv');pv = pv.pv;
[x pid] = max([pv.kVA]);

dt = 30/3600/24;
t = dt:dt:1;
% create one more PV loadshape
idx = length(c.load) + 1;
c.loadshape(idx) = dssloadshape;
nm = 'pv';
c.loadshape(idx).Name = nm;
c.loadshape(idx).Npts = length(t);
c.loadshape(idx).sInterval = 30;
m = zeros(1,length(t));
% find data that we have in PV profile and assign it to m taking into account the time zone for California and the collected data date of Dec 14th 2012 
% special handling needed to avoid rounding error for matlab's time/date
[x, id2] = ismember(round((t+ datenum('2012-12-14 00:00:00'))*24*3600), round((pv(pid).time-8/24)*24*3600) );
m(x) = pv(pid).gi/1000;
% quality control: set NaN values to 0 to avoid fricking out OpenDSS
m(isnan(m)) = 0;
c.loadshape(idx).Mult = m;

% match PV systems in the circuit with calculated profiles
% load new pv system
pv = load('data/scenario2_pv.mat'); 
c.pvsystem = pv.pv;
for i = 1:length(c.pvsystem)
    % assign newly created profile to according PV system
    c.pvsystem(i).daily = nm;
end

%% Save case 1
p = dsswrite(c,'f520',1,[temp '/f520_case1']);
save([temp '/scenario3_case1.mat'],'c','glc');

%% END: Simulation 
Commands={};
SimName = ['Scenario 3 case 1'];
stepsize = '30s';
PVsyst = {'pv'};
SimPer = 24;
NumPts = 2880;
temp = 'tmp_scnr3';
p = dsswrite(c,'f520',1,[temp '/f520_case1']);
DailySimulations(Commands, SimName, stepsize, SimPer, p, PVsyst, NumPts);

%% Case 2: Disaggregate all the loads + a single pv profile. That will be our disaggregated load case. 30s interval
%% loading ...
clear; 
temp = 'tmp_scnr3';
load('data/f520.mat'); load('data/load.mat');

%% Create load profiles 
%deaggregated loadshapes
l1 = ld;
[x idx] = ismember({l1.Name}',strcat('Load.',{c.load.Name}'));

dt = 30/3600/24;
t = dt:dt:1;

offset = 0;%length(c.loadshape);
for i = 1:length(l1)
    % get id of loadshape considering the offset
    id = offset + i;
    c.loadshape(id) = dssloadshape;
    c.loadshape(id).Name = ['loadshape_' c.load(idx(i)).Name];
    c.loadshape(id).sInterval = 30;
    c.loadshape(id).Npts = length(t);
    % set up multiplier; assume the load profile starts at 12am in the morning
    m = [l1(i).x12_00_AM l1(i).x1_00_AM l1(i).x2_00_AM l1(i).x3_00_AM l1(i).x4_00_AM l1(i).x5_00_AM...
        l1(i).x6_00_AM l1(i).x7_00_AM l1(i).x8_00_AM l1(i).x9_00_AM l1(i).x10_00_AM l1(i).x11_00_AM l1(i).x12_00_N...
        l1(i).x1_00_PM l1(i).x2_00_PM l1(i).x3_00_PM l1(i).x4_00_PM l1(i).x5_00_PM l1(i).x6_00_PM...
        l1(i).x7_00_PM l1(i).x8_00_PM l1(i).x9_00_PM l1(i).x10_00_PM l1(i).x11_00_PM l1(i).x12_00_AM];
    % interpolate to get 30 second interval data
    c.loadshape(id).Mult = interp1(0:1/24:1,m,t,'cubic');
    c.load(idx(i)).Daily = c.loadshape(id).Name;
end

%% apply pv profile of the 1MW PV site at hospital to all PV systems
pv = load('Fallbrook_Scenario2/scenario2_pvProfiles.mat','pv');pv = pv.pv;
[x pid] = max([pv.kVA]);

dt = 30/3600/24;
t = dt:dt:1;
% create one more PV loadshape
idx = length(c.load) + 1;
c.loadshape(idx) = dssloadshape;
nm = 'pv';
c.loadshape(idx).Name = nm;
c.loadshape(idx).Npts = length(t);
c.loadshape(idx).sInterval = 30;
m = zeros(1,length(t));
% find data that we have in PV profile and assign it to m taking into account the time zone for California and the collected data date of Dec 14th 2012 
% special handling needed to avoid rounding error for matlab's time/date
[x, id2] = ismember(round((t+ datenum('2012-12-14 00:00:00'))*24*3600), round((pv(pid).time-8/24)*24*3600) );
m(x) = pv(pid).gi/1000;
% quality control: set NaN values to 0 to avoid fricking out OpenDSS
m(isnan(m)) = 0;
c.loadshape(idx).Mult = m;

% match PV systems in the circuit with calculated profiles
% load new pv system
pv = load('data/scenario2_pv.mat'); 
c.pvsystem = pv.pv;
for i = 1:length(c.pvsystem)
    % assign newly created profile to according PV system
    c.pvsystem(i).daily = nm;
end

%% END case 2
p = dsswrite(c,'f520',1,[temp '/f520_case2']);
save([temp '/scenario3_case2.mat'],'c','glc');

%% Case 3: Disaggregate all the pv profiles but keep the loads aggregated (just like case 1). That will be our disaggregated pv case. 30s interval
%% loading ...
clear; 
temp = 'tmp_scnr3';
load('data/f520.mat'); load('data/load.mat');

% load new pv system
pv = load('data/scenario2_pv.mat'); 
c.pvsystem = pv.pv;

%% Create load profiles 
%aggregated loadshapes
l1 = la;
[x idx] = ismember({l1.Name}',strcat('Load.',{c.load.Name}'));

dt = 30/3600/24;
t = dt:dt:1;

offset = 0;%length(c.loadshape);
for i = 1:length(l1)
    % get id of loadshape considering the offset
    id = offset + i;
    c.loadshape(id) = dssloadshape;
    c.loadshape(id).Name = ['loadshape_' c.load(idx(i)).Name];
    c.loadshape(id).sInterval = 30;
    c.loadshape(id).Npts = length(t);
    % set up multiplier; assume the load profile starts at 12am in the morning
    m = [l1(i).x12_00_AM l1(i).x1_00_AM l1(i).x2_00_AM l1(i).x3_00_AM l1(i).x4_00_AM l1(i).x5_00_AM...
        l1(i).x6_00_AM l1(i).x7_00_AM l1(i).x8_00_AM l1(i).x9_00_AM l1(i).x10_00_AM l1(i).x11_00_AM l1(i).x12_00_N...
        l1(i).x1_00_PM l1(i).x2_00_PM l1(i).x3_00_PM l1(i).x4_00_PM l1(i).x5_00_PM l1(i).x6_00_PM...
        l1(i).x7_00_PM l1(i).x8_00_PM l1(i).x9_00_PM l1(i).x10_00_PM l1(i).x11_00_PM l1(i).x12_00_AM];
    % interpolate to get 30 second interval data
    c.loadshape(id).Mult = interp1(0:1/24:1,m,t,'cubic');
    c.load(idx(i)).Daily = c.loadshape(id).Name;
end

%% deaggregate PV profiles
pv = load('Fallbrook_Scenario2/scenario2_pvProfiles.mat','pv');pv = pv.pv;
% match PV systems in the circuit with calculated profiles
[x id] = ismember({c.pvsystem.Name}',{pv.Name}');
dt = 30/3600/24;
t = dt:dt:1;
offset = length(c.load);
for i = 1:length(id)
    c.loadshape(offset+i) = dssloadshape;
    % loadshape name
    nm = ['pv' sprintf('%d',i)];
    c.loadshape(offset+i).Name = nm;
    c.loadshape(offset+i).Npts = length(t);
    c.loadshape(offset+i).sInterval = 30;
    m = zeros(1,length(t));
    % find data that we have in PV profile and assign it to m taking into account the time zone for California and the collected data date of Dec 14th 2012 
    % special handling needed to avoid rounding error for matlab's time/date
    [x, id2] = ismember(round((t+ datenum('2012-12-14 00:00:00'))*24*3600), round((pv(id(i)).time-8/24)*24*3600) );
    m(x) = pv(id(i)).gi/1000;
    % quality control: set NaN values to 0 to avoid fricking out OpenDSS
    m(isnan(m)) = 0;
    c.loadshape(offset+i).Mult = m;
    
    % assign newly created profile to according PV system
    c.pvsystem(i).daily = nm;
end

%% END case 3
p = dsswrite(c,'f520',1,[temp '/f520_case3']);
save([temp '/scenario3_case3.mat'],'c','glc');

%% Case 4: Disaggregate both (i.e. loads and pv). That will be our disaggregated loads and PVs case. 30s interval
%% loading ...
clear; 
temp = 'tmp_scnr3';
load('data/f520.mat'); load('data/load.mat');
% load new pv system
pv = load('data/scenario2_pv.mat'); 
c.pvsystem = pv.pv;

%% Create load profiles 
%deaggregated loadshapes
l1 = ld;
[x idx] = ismember({l1.Name}',strcat('Load.',{c.load.Name}'));

dt = 30/3600/24;
t = dt:dt:1;

offset = 0;%length(c.loadshape);
for i = 1:length(l1)
    % get id of loadshape considering the offset
    id = offset + i;
    c.loadshape(id) = dssloadshape;
    c.loadshape(id).Name = ['loadshape_' c.load(idx(i)).Name];
    c.loadshape(id).sInterval = 30;
    c.loadshape(id).Npts = length(t);
    % set up multiplier; assume the load profile starts at 12am in the morning
    m = [l1(i).x12_00_AM l1(i).x1_00_AM l1(i).x2_00_AM l1(i).x3_00_AM l1(i).x4_00_AM l1(i).x5_00_AM...
        l1(i).x6_00_AM l1(i).x7_00_AM l1(i).x8_00_AM l1(i).x9_00_AM l1(i).x10_00_AM l1(i).x11_00_AM l1(i).x12_00_N...
        l1(i).x1_00_PM l1(i).x2_00_PM l1(i).x3_00_PM l1(i).x4_00_PM l1(i).x5_00_PM l1(i).x6_00_PM...
        l1(i).x7_00_PM l1(i).x8_00_PM l1(i).x9_00_PM l1(i).x10_00_PM l1(i).x11_00_PM l1(i).x12_00_AM];
    % interpolate to get 30 second interval data
    c.loadshape(id).Mult = interp1(0:1/24:1,m,t,'cubic');
    c.load(idx(i)).Daily = c.loadshape(id).Name;
end

%% deaggregate PV profiles
pv = load('Fallbrook_Scenario2/scenario2_pvProfiles.mat','pv');pv = pv.pv;
% match PV systems in the circuit with calculated profiles
[x id] = ismember({c.pvsystem.Name}',{pv.Name}');
dt = 30/3600/24;
t = dt:dt:1;
offset = length(c.load);
for i = 1:length(id)
    c.loadshape(offset+i) = dssloadshape;
    % loadshape name
    nm = ['pv' sprintf('%d',i)];
    c.loadshape(offset+i).Name = nm;
    c.loadshape(offset+i).Npts = length(t);
    c.loadshape(offset+i).sInterval = 30;
    m = zeros(1,length(t));
    % find data that we have in PV profile and assign it to m taking into account the time zone for California and the collected data date of Dec 14th 2012 
    % special handling needed to avoid rounding error for matlab's time/date
    [x, id2] = ismember(round((t+ datenum('2012-12-14 00:00:00'))*24*3600), round((pv(id(i)).time-8/24)*24*3600) );
    m(x) = pv(id(i)).gi/1000;
    % quality control: set NaN values to 0 to avoid fricking out OpenDSS
    m(isnan(m)) = 0;
    c.loadshape(offset+i).Mult = m;
    
    % assign newly created profile to according PV system
    c.pvsystem(i).daily = nm;
end

%% END case 4
p = dsswrite(c,'f520',1,[temp '/f520_case4']);
save([temp '/scenario3_case4.mat'],'c','glc');

%% Case 5: Disaggregate both but instead of using 30 seconds data use 1h data. That will be our low resolution case. After that we run OpenDSS and save the results.
%% loading ...
clear; 
temp = 'tmp_scnr3';
load('data/f520.mat'); load('data/load.mat');
% load new pv system
pv = load('data/scenario2_pv.mat'); 
c.pvsystem = pv.pv;

%% Create load profiles 
%deaggregated loadshapes
l1 = ld;
[x idx] = ismember({l1.Name}',strcat('Load.',{c.load.Name}'));

dt = 1/24;
t = dt:dt:1;

offset = 0;%length(c.loadshape);
for i = 1:length(l1)
    % get id of loadshape considering the offset
    id = offset + i;
    c.loadshape(id) = dssloadshape;
    c.loadshape(id).Name = ['loadshape_' c.load(idx(i)).Name];
    c.loadshape(id).interval = 1;
    c.loadshape(id).Npts = length(t);
    % set up multiplier; assume the load profile starts at 1am in the morning
    m = [l1(i).x1_00_AM l1(i).x2_00_AM l1(i).x3_00_AM l1(i).x4_00_AM l1(i).x5_00_AM...
        l1(i).x6_00_AM l1(i).x7_00_AM l1(i).x8_00_AM l1(i).x9_00_AM l1(i).x10_00_AM l1(i).x11_00_AM l1(i).x12_00_N...
        l1(i).x1_00_PM l1(i).x2_00_PM l1(i).x3_00_PM l1(i).x4_00_PM l1(i).x5_00_PM l1(i).x6_00_PM...
        l1(i).x7_00_PM l1(i).x8_00_PM l1(i).x9_00_PM l1(i).x10_00_PM l1(i).x11_00_PM l1(i).x12_00_AM];
    % interpolate to get 30 second interval data
    c.loadshape(id).Mult = m;
    c.load(idx(i)).Daily = c.loadshape(id).Name;
end

%% deaggregate pv profiles
pv = load('Fallbrook_Scenario2/scenario2_pvProfiles.mat','pv');pv = pv.pv;
% match PV systems in the circuit with calculated profiles
[x id] = ismember({c.pvsystem.Name}',{pv.Name}');
dt = 1/24;
t = dt:dt:1;
offset = length(c.load);
for i = 1:length(id)
    c.loadshape(offset+i) = dssloadshape;
    % loadshape name
    nm = ['pv' sprintf('%d',i)];
    c.loadshape(offset+i).Name = nm;
    c.loadshape(offset+i).Npts = length(t);
    c.loadshape(offset+i).interval = 1;
    m = zeros(1,length(t));
    % find data that we have in PV profile and assign it to m taking into account the time zone for California and the collected data date of Dec 14th 2012 
    % special handling needed to avoid rounding error for matlab's time/date
    [x, id2] = ismember(round((t+ datenum('2012-12-14 00:00:00'))*24), round((pv(id(i)).time-8/24)*24) );
    id2 = id2(id2>0);
    m(x) = pv(id(i)).gi(id2)/1000;
    % quality control: set NaN values to 0 to avoid fricking out OpenDSS
    m(isnan(m)) = 0;
    c.loadshape(offset+i).Mult = m;
    
    % assign newly created profile to according PV system
    c.pvsystem(i).daily = nm;
end

%% END case 5
p = dsswrite(c,'f520',1,[temp '/f520_case5']);
save([temp '/scenario3_case5.mat'],'c','glc');

%% Case 6: Disaggregate all the loads + a single pv profile on a clear San Diego, Dec 14th 2013, day. 30s interval
%% loading ...
clear; 
temp = 'tmp_scnr3';
load('data/f520.mat'); load('data/load.mat');
% load clear sky irradiance
load('data/20121214_clrSkyGHI.mat');

%% Create load profiles 
%deaggregated loadshapes
l1 = ld;
[x idx] = ismember({l1.Name}',strcat('Load.',{c.load.Name}'));

dt = 30/3600/24;
t = dt:dt:1;

offset = 0;%length(c.loadshape);
for i = 1:length(l1)
    % get id of loadshape considering the offset
    id = offset + i;
    c.loadshape(id) = dssloadshape;
    c.loadshape(id).Name = ['loadshape_' c.load(idx(i)).Name];
    c.loadshape(id).sInterval = 30;
    c.loadshape(id).Npts = length(t);
    % set up multiplier; assume the load profile starts at 12am in the morning
    m = [l1(i).x12_00_AM l1(i).x1_00_AM l1(i).x2_00_AM l1(i).x3_00_AM l1(i).x4_00_AM l1(i).x5_00_AM...
        l1(i).x6_00_AM l1(i).x7_00_AM l1(i).x8_00_AM l1(i).x9_00_AM l1(i).x10_00_AM l1(i).x11_00_AM l1(i).x12_00_N...
        l1(i).x1_00_PM l1(i).x2_00_PM l1(i).x3_00_PM l1(i).x4_00_PM l1(i).x5_00_PM l1(i).x6_00_PM...
        l1(i).x7_00_PM l1(i).x8_00_PM l1(i).x9_00_PM l1(i).x10_00_PM l1(i).x11_00_PM l1(i).x12_00_AM];
    % interpolate to get 30 second interval data
    c.loadshape(id).Mult = interp1(0:1/24:1,m,t,'cubic');
    c.load(idx(i)).Daily = c.loadshape(id).Name;
end

%% apply pv profile of the 1MW PV site at hospital to all PV systems
pv = load('Fallbrook_Scenario2/scenario2_pvProfiles.mat','pv');pv = pv.pv;
[x pid] = max([pv.kVA]);

dt = 30/3600/24;
t = dt:dt:1;
% create one more PV loadshape
idx = length(c.load) + 1;
c.loadshape(idx) = dssloadshape;
nm = 'pv';
c.loadshape(idx).Name = nm;
c.loadshape(idx).Npts = length(t);
c.loadshape(idx).sInterval = 30;
m = zeros(1,length(t));
% find data that we have in PV profile and assign it to m taking into account the time zone for California and the collected data date of Dec 14th 2012 
% special handling needed to avoid rounding error for matlab's time/date
[x, id2] = ismember(round((t+ datenum('2012-12-14 00:00:00'))*24*3600), round((ghi.time)*24*3600) );
m(x) = ghi.ghi(id2)/1000;
% quality control: set NaN values to 0 to avoid fricking out OpenDSS
m(isnan(m)) = 0;
c.loadshape(idx).Mult = m;

% match PV systems in the circuit with calculated profiles
pv = load('data/scenario2_pv.mat'); 
c.pvsystem = pv.pv;
for i = 1:length(c.pvsystem)
    % assign newly created profile to according PV system
    c.pvsystem(i).daily = nm;
end

%% END case 6
p = dsswrite(c,'f520',1,[temp '/f520_case6']);
save([temp '/scenario3_case6.mat'],'c','glc');