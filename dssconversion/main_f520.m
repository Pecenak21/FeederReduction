%% summary
% use dssconversion to create initial feeder
% modify and validate circuit
% output: data/f520.mat
%% Main file for processing data with Feeder 520 
% Get data
% if you want to specify local directory then use following codes (try not to commit your local path, you can browse to the directory and use 'pwd' instead):
% strDirectory='D:\Sync Software\Other Code\SystemConversion\Convert2dss';
% strDirectory='C:\Work\Projects\2012\1787-UCSD PV\Simulation\System\SDGE';
% cd(strDirectory); 
curdir = pwd;
fn = [curdir '/dssconversion/custdata/520ForENERNEX.xlsx'];
d = excel2obj(fn);
% given linecode
glc = excel2obj( 'linecode.xlsx' );
glc = glc.LineCode;

% convert to dss struct
c = dssconversion( d, glc ); c_bk = c; 

% Change Circuit settings
c.circuit.bus1 = 'sourcebus';
% c.circuit.pu = 1.00;% changed from 1.00
c.circuit.basekv = 12.47;
c.basevoltages = [12.47 12];
% substation
c.transformer(end+1) = c.transformer(end);
c.transformer(end).Name = 'AVOCADO';
c.transformer(end).Buses = {'SourceBus' '0520'};
c.transformer(end).Conns = {'delta' 'wye'};
c.transformer(end).kVs = [12.47 12];
c.transformer(end).kVAs = [28000 28000];% changed from [10000 10000]
c.transformer(end).XHL = 1.0871;
c.transformer(end).sub = 'y';
c.transformer(end).Rs = [0.103 0.103];
c.regcontrol(end+1) = c.regcontrol(end);
c.regcontrol(end).Name = 'Avocado';
c.regcontrol(end).transformer = 'Avocado';
c.regcontrol(end).vreg = 121;

% change load setting
loads_1Phase = c.load([c.load.Phases]==1);
[a1, b1] = ismember({loads_1Phase.Name}, {c.load.Name});
for i=1:length(b1)
    c.load(b1(i)).Kv = 6.9282;
    c.load(b1(i)).Kvar = c.load(b1(i)).Kvar*2.5;
end
loads_3Phase = c.load([c.load.Phases]==3);
[a3, b3] = ismember({loads_3Phase.Name}, {c.load.Name});
for j=1:length(b3)
    c.load(b3(j)).Kv = 12;
    if c.load(b3(j)).Kvar < 189 % the 4 biggest loads stay unchanged
       c.load(b3(j)).Kvar = c.load(b3(j)).Kvar*2.5;
    end
end

% Optimizing capacitor banks to match reactive power 
% Modify existing caps
c.capacitor(1).Kvar = 1420;% 2900/sqrt(3);
c.capacitor(1).Numsteps = 20;
c.capacitor(1).Name = ['cap_' c.capacitor(1).Name];
c.capcontrol(1).Capacitor = c.capacitor(1).Name;
c.capacitor(2).Kvar = 1400;%650/sqrt(3);
c.capacitor(2).Numsteps = 20;
c.capacitor(2).Name = ['cap_' c.capacitor(2).Name];
c.capcontrol(2).Capacitor = c.capacitor(2).Name;

% Adding cap banks at appropriate locations
c.capacitor(3) = c.capacitor(2);
c.capacitor(3).Bus1 = '05201643A.1.2.3';
c.capacitor(3).Name = ['cap_' cleanBus(c.capacitor(3).Bus1) 'var'];
c.capacitor(3).kvar = 1300;
c.capacitor(3).Numsteps = 25;

c.capacitor(4) = c.capacitor(2);
c.capacitor(4).Bus1 = '05201947';
c.capacitor(4).Name = ['cap_' c.capacitor(4).Bus1 'var'];
c.capacitor(4).kvar = 120;

c.capacitor(5) = c.capacitor(2);
c.capacitor(5).Bus1 = '05201349';
c.capacitor(5).Name = ['cap_' c.capacitor(5).Bus1 'var'];
c.capacitor(5).kvar = 60;

% %% Adding/Editing CapControls 
c.capcontrol(1).Capacitor = c.capacitor(1).Name;
c.capcontrol(1).Vmax = 126;

c.capcontrol(2).Capacitor = c.capacitor(2).Name;
c.capcontrol(2).Vmax = 126;

% %% Changing Generator Settings
c.generator.Kv = 12;
c.generator.Model = 3;

% %% Changing Transformer Settings
for i=1:length(c.transformer)
    c.transformer(i).Conns = {'wye' 'wye'};
    c.transformer(i).Windings = 2;
    %was [5000 5000]-made it bigger for less losses. the reactive power 
    %validation with SynerGEE does no match anymore (when set to 20k). 
    %the reason being that the cap banks are not big enough to support bigger trafos. However,
    %for the daily runs this size shows better voltage levels, while the
    %cap banks don't do anything anyway. 
    c.transformer(i).kVAs = [20000 20000];
    c.regcontrol(i).vreg = 122;
    c.regcontrol(i).band = 1;
    c.regcontrol(i).ptratio = 60;
    c.regcontrol(i).EventLog = 'yes';
    c.regcontrol(i).PTPhase = [];
    c.regcontrol(i).vlimit = [];
    c.regcontrol(i).revNeutral = [];
    c.transformer(i).Taps = [];
end

%push the voltage regulator to a different bus. it is important for daily
%simulations. at the end of the feeder the voltage at bus B drops below 93%
%around 8pm. by pushing the regulator back the voltage in this region can
%be boosted up to stay above 95%
c.transformer(3).Buses={'05201400' '05201401'};
[x, y] = ismember('05201400_05201401', c.line.Name(:));
n = c.line(y);
c.line(y) = [];%remove the line to make room for the transformer
n.bus1 = '05201438.1.2.3';
n.bus2 = '05201438A.1.2.3';
n.Name = '05201438_05201438A';
c.line(end+1) = n;%place a line a the former location of the regulator

% %% add "sourcebus" node to buslist
c.buslist.id = ['SourceBus'; c.buslist.id];
c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];

% %% addPV to the circuit
% this parts is not yet optimized - loading pvsystems from an existing
% manually generated profile.
pv = dssparse(which('f520_pvsystem_pvsystems.dss'));
c.pvsystem = pv.pvsystem;

% %% Close all switches
% seems to be a line to 'nowhere'. when open - voltage is 0V and screws the voltage plots. Does not change
% the overall results in either mode.
c.switch(13).Action = 'Close';

% %% Daily Simulations - modification of the original circuit
% Replace the Generator with pvSystems (the generator is useful only for the power
% validation task since that is how it is apparently done in the SynerGee simulation)
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

% %% Changing Regcontrol settings to better suit the daily simulaitons. The
% previous settings were necessary to validate active and reactive powers.
for i=1:length(c.regcontrol)
    c.regcontrol(i).vreg = 121;
    c.regcontrol(i).band = 1;
    c.regcontrol(i).ptratio = 57.75;
    c.regcontrol(i).EventLog = 'yes';
    c.regcontrol(i).PTPhase = 3;
    c.regcontrol(i).vlimit = 126;
    c.regcontrol(i).revNeutral = [];
end
c.regcontrol(end).vreg = 122;%higher voltage at substation
c.regcontrol(end).band = 2;%higher control bandwith at substation

% %% Changing Capcontrol settings to better suit the daily simulaitons.
for i=1:length(c.capcontrol)
    c.capcontrol(i).OFFsetting = 125;
    c.capcontrol(i).ONsetting = 119;
    c.capcontrol(i).Vmin = 117;
    c.capcontrol(i).VoltOverride = 'TRUE';
    c.capcontrol(i).EventLog = 'yes';
end
%%
p=dsswrite(c);
c.load(1731).Kw = 525;
c.load(1732).Kw = 1313;
validatepower(3,'kvar',0,'',[1 0 0 0 0],c,glc,0);
validatepower(3,'kw',0,'',[1 0 0 0 0],c,glc,0);
faultstudy(c,d,p);
% %% Save original set up
save('tmp/f520.mat','c','p','glc');