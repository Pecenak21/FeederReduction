function c = circuit_changes(c)
% c = circuit_changes(c)
%
% PURPOSE : Tweak the circuit that has been converted to OpenDSS format.
% This is sometimes necessary (1) to match load flow and short circuit
% results and (2) to add devices/information that is not included in the
% original system.
%
%
% INPUT :   c: unmodified circuit in OpenDSS format
%
% OUTPUT :  c: modified circuit in OpenDSS format


%% Optimizing capacitor banks (Be careful when combining cases!) 
%NOTE: It has been tested and the best result (matching with the given power flow data from utility) is to combine cases 1 and 2. 
% Case 1. use modified caps
c.capacitor(1).kvar = 1600;% 2900/sqrt(3);
c.capacitor(2).kvar = 1500;%650/sqrt(3);

% Case 2. create new caps
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
for i=1:1:length(b1)
    c.load(b1(i)).Kv = 6.9282;
    c.load(b1(i)).Kvar = c.load(b1(i)).Kvar*2.5;
end
loads_3Phase = c.load([c.load.Phases]==3);
[a3 b3] = ismember({loads_3Phase.Name}, {c.load.Name});
for j=1:1:length(b3)
    c.load(b3(j)).Kv = 12;
    if c.load(b3(j)).Kvar < 189 % the 4 biggest loads stay unchanged
       c.load(b3(j)).Kvar = c.load(b3(j)).Kvar*2.5;
    end
end

%% Changing Transformer Settings
for k=1:1:length(c.transformer)
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
if isfield(c,'switch')
    c.switch(13).Action = 'Close';
end