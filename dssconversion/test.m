clear all
clc

%% Main file for processing data with Feeder 520 
% Get data. If your working folder is not in this repo, change it accordingly on below line
% strDirectory='C:\Work\Projects\2012\1787-UCSD PV\Simulation\System\SDGE';
% cd(strDirectory); 
curdir = pwd;
fn = [curdir '/custdata/520ForENERNEX.xlsx'];
d = excel2obj(fn);
% given linecode
glc = excel2obj( 'linecode.xlsx' );
glc = glc.LineCode;

% convert to dss struct
c = dssconversion( d, glc );

%% Optimizing capacitor banks to match reactive power
% Modify existing caps
c.capacitor(1).kvar = 1600;% 2900/sqrt(3);
c.capacitor(2).kvar = 1500;%650/sqrt(3);

% Adding cap banks at appropriate locations
c.capacitor(3) = c.capacitor(2);
c.capacitor(3).Bus1 = '05201643A';
c.capacitor(3).Name = ['cap_' c.capacitor(3).Bus1];
c.capacitor(3).kvar = 1230;

c.capacitor(4) = c.capacitor(3);
c.capacitor(4).Bus1 = '05201644A';
c.capacitor(4).Name = ['cap_' c.capacitor(4).Bus1];
c.capacitor(4).kvar = 100;

%% Changing Generator Settings
c.generator.Kv = 12;
c.generator.Model = 3;

%% Changing Load Settings
loads_1Phase = c.load([c.load.Phases]==1);
[a1 b1] = ismember({loads_1Phase.Name}, {c.load.Name});
for i=1:length(b1)
    c.load(b1(i)).Kv = 6.9282;
    c.load(b1(i)).Kvar = c.load(b1(i)).Kvar*2.5;
end
loads_3Phase = c.load([c.load.Phases]==3);
[a3 b3] = ismember({loads_3Phase.Name}, {c.load.Name});
for j=1:length(b3)
    c.load(b3(j)).Kv = 12;
    if c.load(b3(j)).Kvar < 189 % the 4 biggest loads stay unchanged
       c.load(b3(j)).Kvar = c.load(b3(j)).Kvar*2.5;
    end
end

%% Changing Transformer Settings
for k=1:length(c.transformer)
    c.transformer(k).Conns = {'wye' 'wye'};
    c.transformer(k).Windings = 2;
    c.transformer(k).NumTaps = [];
    c.transformer(k).kVAs = [5000 5000];
    c.regcontrol(k).vreg = 125;
    c.regcontrol(k).band = 3;
    c.regcontrol(k).ptratio = 60;
    c.regcontrol(k).EventLog = 'yes';
    c.regcontrol(k).PTPhase = [];
    c.regcontrol(k).vlimit = [];
    c.regcontrol(k).revNeutral = [];
end
% substation
c.transformer(end+1) = c.transformer(end);
c.transformer(end).Name = 'AVOCADO';
c.transformer(end).Buses = {'SourceBus' '0520.1.2.3.0'};
c.transformer(end).Conns = {'delta' 'wye'};
c.transformer(end).kVs = [115 12];
c.transformer(end).kVAs = [28000 28000];% changed from [10000 10000]
c.transformer(end).XHL = 1.0871;
c.transformer(end).sub = 'y';
c.transformer(end).Rs = [0.103 0.103];
c.regcontrol(end+1) = c.regcontrol(end);
c.regcontrol(end).Name = 'Avocado';
c.regcontrol(end).transformer = 'Avocado';
c.regcontrol(end).vreg = 121;

%% addPV to the circuit
% this parts is not yet optimized - loading pvsystems from an existing
% manually generated profile.
pv = dssparse('data/f520_pvsystem_pvsystems.dss');
c.pvsystem = pv.pvsystem;
% c = addPV(c,'custdata/520_PV.xlsx');
% l=1;
% while l <= length(c.pvsystem)
%     if c.pvsystem(l).kVA > 900 %remove hospital pv systems
%         c.pvsystem(l) =[];continue;end
%     c.pvsystem(l).kv = 12;
%     c.pvsystem(l).pf = 1;
%     c.pvsystem(l).phases = 3;
%     c.pvsystem(l).conn = 'wye';
%     l=l+1;
% end

%% Change Circuit settings
c.circuit.bus1 = '';
c.circuit.pu = 1.00;% changed from 1.05
c.circuit.basekv = [];
c.basevoltages = [115 12];

%% Close all switches
% seems to be a line to 'nowhere'. when open - voltage is 0V and screws the voltage plots. Does not change
% the overall results in either mode.
c.switch(13).Action = 'Close';

%% Replace the Generator with pvSystems (the generator is useful only for the power validation task since that is how it is apparently done in the SynerGee simulation)
c.generator.Enabled = 'no';
PV = c.pvsystem(end);
PV.Name = 'CG199999_40477_Calle_Roxanne';
PV.bus1 = '05201644A';
PV.phases = 3;
PV.kv = 12;
PV.kVA = 999;
c.pvsystem(end+1) = PV;
PV.Name = 'CG199999_40177_Calle_Roxanne';
PV.bus1 = '05201643A';
PV.kVA = 998.9;
c.pvsystem(end+1) = PV;

%% Save circuit
if ~exist('tmp','dir'),  mkdir('tmp'); end;
p = dsswrite(c,'f520',1,'tmp/f520');
save('tmp/f520.mat');

%% Run simulation to test matching with given data
% run if kw or kvar validation desired
validatepower(3,'kw',0,'',[1 0 0 0 0],c,glc,0,p);
validatepower(3,'kvar',0,'',[1 0 0 0 0],c,glc,0,p);

%% Next time just load data
% load tmp/f520.mat;

%% Daily Simulations/ Quasi-state solution for a day; Create loadshapes for loads and pvprofiles
%% Get load profiles
ls = dssparse('data/Loadshapes.dss');
c.loadshape = ls.loadshape;

%% plot to check load shapes
figure, hold on
cmap = jet(256);
for i = 1:length(c.loadshape)
    plot([c.loadshape(i).Mult],'color',cmap(round(i/length(c.loadshape)*255),:,:),'linewidth',2);
end
xlim([1 24]);
xlabel('Time of Day, [hour]','fontsize',20);
ylabel('Mulitplier, [1]','fontsize',20);
box on, grid on; 
set(gcf,'color','w','position',[50 50 1200 900]); set(gca,'fontsize',20);

% NOTE: Profile 9 is a small commercial or small industrial load with the load consistently at around .5 power unit. We could use it for the hospital. 
%% time set up for 24 hours of a day in 30 second interval on Dec 14th
dt = 30/3600/24;
t = 0:dt:24-dt;

%% create load shape from given ones; interpolate to get 30s interval data
for i = 1:length(c.loadshape)
    c.loadshape(i).Interval = [];
    c.loadshape(i).sInterval = 30;
    c.loadshape(i).Npts = length(t);
    % set up multiplier
    m = ls.loadshape(i).Mult;
    % assume the load profile starts at 12am in the morning
    m = [m m(1)];
    % interpolate to get 30 second interval data
    m = interp1(0:24,m,t,'cubic');
    c.loadshape(i).Mult = m;
end

%% Assign a random load profile out of 15 given loadshapes (exclude loadshape #9 for comercial/industrial loads) to each load
for m=1:length(c.load)
    % get a random load out of the 14 ligit ones
    id = ceil(rand(1,1)*14);
    % skip loadshape #9 by adding 1 to the index
    if id >= 9, id = id + 1; end
    c.load(m).Daily = c.loadshape(id).Name;
    c.load(m).Status = 'variable';
    c.load(m).Conn = 'wye';
end
%% assign loadshape #9 to the hospital load
% search for the hospital loads. There are 2 loads at the hospital > 500 kW
[x, id] = find([c.load.Kw]' > 500);
% assign loadshape #9 to these loads
for i = x
    c.load(m).Daily = c.loadshape(9).Name;
end

%% Load PV deaggregated profiles and generate PV profiles accordingly
pv = load('PVProfilesFeeder520/pv.mat','pv');pv = pv.pv;
% match PV systems in the circuit with calculated profiles
[x id] = ismember({c.pvsystem.Name}',{pv.Name}');
for i = 1:length(id)
    c.loadshape(end+1) = dssloadshape;
    % loadshape name
    nm = ['pv' sprintf('%d',i)];
    c.loadshape(end).Name = nm;
    c.loadshape(end).Npts = length(t);
    c.loadshape(end).sInterval = 30;
    m = zeros(1,length(t));
    % find data that we have in PV profile and assign it to m taking into account the time zone for California and the collected data date of Dec 14th 2012 
    [x, id2] = ismember(t+ datenum([2012 12 14 00 00 00]), pv(id(i)).time-8/24 );
    m(x) = pv(id(i)).gi(id2)/1000;
    c.loadshape(end).Mult = m;
    
    % assign newly created profile to according PV system
    c.pvsystem(i).daily = nm;
end

%% Add EnergyMeter at Substation
n = dssenergymeter('Name', 'Sub');
n.element = 'Line.0520_05201';
n.terminal = 1;
c.energymeter = n;

%% Place a monitor at the Substation low-site
n = dssmonitor('Name','Sub');
n.Element = 'Line.0520_05201';
n.Mode = 0;
c.monitor = n;

%% save edited circuit
p = dsswrite(c,'f520',1,'tmp/f520');

%% Run Simulations

Commands = {};%any additional changes to the circuit on the fly (e.g. settign values, etc.)
Commands = ['Loadshape.Dom1.Mult = (' num2str(ones(1,24)) ')'];
Commands = {Commands; ['Loadshape.PVProfile.Mult = (' num2str(ones(1,24)) ')']};
% Commands = ['Loadshape.Dom1.Mult = (' num2str(ones(1,24)) ')'];
% Commands = {Commands; ['Loadshape.PVProfile.Mult = (' num2str(zeros(1,24)) ')']};
% Commands = ['Loadshape.Dom1.Mult = (' num2str(ShapeLoad) ')'];
% Commands = {Commands; ['Loadshape.PVProfile.Mult = (' num2str(ShapePV) ')']};

SimName = 'Load 100% & PV 100%';%simulation title, will appear on the figures and figure names when saved
stepsize = '1h';%time will advance based on stepsize (e.g. 1h or 15m) for each solve
NumPts = 24;%number of 'solves' (e.g. 24 for a daily simulaiton with 1h stepsize)

% 3D plots parameters:
z_saturation=[0.70 1.15];% z_saturation influences color spectrum, the higher the range, the farther away from extreme color
z_range=(0.75:0.05:1.05);% range of values displayed, also determines z-axis numbering

DailySimulations(Commands, SimName, stepsize, NumPts, p, PVSys, z_saturation, z_range);


