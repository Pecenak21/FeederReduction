function [c1] = circuit_changes_Durox(c1, node_bk)

% line impedances need to have tiny values rather than 0s
c1.linecode(87).R1 = 0.000001;
c1.linecode(87).X1 = 0.000001;
c1.linecode(87).R0 = 0.000001;
c1.linecode(87).X0 = 0.000001;

%% change impedance of the line codes (it appears that the impedances are in pu rather than in ohms)
Vbase = 12470;
Mbase = 100000000;
Zbase = Vbase^2/Mbase;
for i=1:length(c1.linecode)
    c1.linecode(i).R1 = c1.linecode(i).R1*(Zbase);
    c1.linecode(i).X1 = c1.linecode(i).X1*(Zbase);
    c1.linecode(i).R0 = c1.linecode(i).R0*(Zbase);
    c1.linecode(i).X0 = c1.linecode(i).X0*(Zbase);
    c1.linecode(i).C1 = c1.linecode(i).C1*(Zbase);
    c1.linecode(i).C0 = c1.linecode(i).C0*(Zbase);
end

%% Add EnergyMeter at Substation
n = dssenergymeter('Name', 'Sub');
n.element = 'Line.05465_2308871E_05465_1_05465';
n.terminal = 1;
c1.energymeter = n;

%% add remaining PV manually
PV = dsspvsystem('Name', '44');
PV.bus1 = '34';
PV.phases = 3;
PV.kv = 0.208;
PV.kVA = 285;
PV.pf = 1;
PV.EffCurve = 'eff';
PV.PTCurve = 'PvsT';
PV.Tdaily = 'Temp';
PV.Pmpp = 285;
PV.cutin = 0.01;
PV.cutout = 1.2;
PV.conn = 'delta';
PV.irradiance = 1;
PV.temperature = 25;
c1.pvsystem = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '45';
PV.bus1 = '38';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '46';
PV.bus1 = '39';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '47';
PV.bus1 = '40';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '48';
PV.bus1 = '41';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '49';
PV.bus1 = '42';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '53';
PV.bus1 = '43';
c1.pvsystem(end+1) = PV;

%% add efficiency, temp, and PvsT curves
eff = dssxycurve('Name', 'eff');
eff.npts = 4;
eff.Xarray = [0.01	0.2	0.4	1.0];
eff.Yarray = [1.0	1	1	1.0];
c1.xycurve = eff;

PvsT = dssxycurve('Name', 'PvsT');
PvsT.npts = 4;
PvsT.Xarray = [0	25	75	100];
PvsT.Yarray = [1.0	1	1	1.0];
c1.xycurve(end+1) = PvsT;

temp = dsstshape('Name', 'Temp');
temp.npts = 24;
temp.temp = [25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25];
c1.tshape = temp;

%% remove additional lines (during conversion switches have been added that need to be replaced with PVs)
PVs=c1.pvsystem(:).Name;
[aa,bb]=ismember(PVs, c1.line(:).Name);
bb=bb(aa);
c1.line(bb)='';

%% remove additional lines (during conversion all loads were connected to the circuit through switches. These
%% are not required.
[aa,bb] = ismember(cell2mat(c1.line.Phases(:)),0);
Line_bus1 = regexprep(c1.line(aa).bus1,'(\.\d+)+$','');
Line_bus2 = regexprep(c1.line(aa).bus2,'(\.\d+)+$','');
c1.line(aa) = '';% remove the lines
Load_bus1 = regexprep(c1.load.bus1(:),'(\.\d+)+$','');
[aa,bb]=ismember(Load_bus1,Line_bus2);
for i=1:length(c1.load.bus1(:))% change loads' bus1 from the end of the line to the beginning of it.
    ind = bb(i);
    ph = c1.load(i).Name(end);
    if strcmpi(ph,'X')
        c1.load(i).bus1 = [Line_bus1{ind} '.1'];
    elseif strcmpi(ph,'Y')
        c1.load(i).bus1 = [Line_bus1{ind} '.2'];
    elseif strcmpi(ph,'Z')
        c1.load(i).bus1 = [Line_bus1{ind} '.3'];
    else
        c1.load(i).bus1 = Line_bus1{ind};
    end
end

%% add missing lines (after conversion three lines are missing - not sure if they are even needed)
n = c1.line(1);
n.Name = 'ND46615856_05465_860469E_05465';
n.LineCode = 'CLP_2C_#2_5/12/16KV';
n.Phases = 2;
n.bus1 = 'ND46615856_05465.1.2';
n.bus2 = '860469E_05465.1.2';
n.Length = 74.371198;
c1.line(end+1) = n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 'ND46616354_05465_1865426E_05465';
n.LineCode = 'CLP_2C_#2_5/12/16KV';
n.Phases = 2;
n.bus1 = 'ND46616354_05465.2.3';
n.bus2 = '1865426E_05465.2.3';
n.Length = 54.863998;
c1.line(end+1) = n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 'ND63836859_05465_2210893E_05465_';
n.LineCode = 'CLP_2C_#2_5/12/16KV';
n.Phases = 2;
n.bus1 = 'ND63836859_05465.1.2';
n.bus2 = '2210893E_05465_1.1.2';
n.Length = 11.5824;
c1.line(end+1) = n;

%% add voltage setting to loads
for i=1:length(c1.load)
    if (c1.load(i).Phases > 1)
        c1.load(i).Kv = 12;
    else
        c1.load(i).Kv = 12/sqrt(3);
    end
    c1.load(i).NumCust = 1;
end

%% remove swiches connecting capacitors. These swiches are from Network.txt and represent sections rather
% than lines. 
buses = c1.capacitor(:).Bus1;
LineBus1 = regexprep(c1.line(:).bus1,'(\.\d+)+$','');
LineBus2 = regexprep(c1.line(:).bus2,'(\.\d+)+$','');
[a,b] = ismember(buses,LineBus2);
CapBus = LineBus1(b);
[aa,bb] = ismember(CapBus,LineBus2);
LineNames = c1.line(bb).Name;% monitor buses upstream from Caps
for i=1:length(b)
    c1.capacitor(i).Bus1 = CapBus(i);
    c1.capcontrol(i).Element = ['line.' LineNames{i}];
end
c1.line(b) = '';

%% change capacitor settings
% c1.capacitor(2).Kv = 12;
c1.capacitor(2).Kv = 12*0.95;
c1.capacitor(2).Name = ['cap_' c1.capacitor(2).Name];
% c1.capacitor(3).Kv = 12;
c1.capacitor(3).Kv = 12*0.95;
c1.capacitor(3).Name = ['cap_' c1.capacitor(3).Name];
c1.capacitor(1) = [];% remove inactive capacitor

%% change capcontrol settings
c1.capcontrol(2).Type = 'voltage';
c1.capcontrol(2).CTRatio = 100;
c1.capcontrol(2).Eventlog = 'True';
c1.capcontrol(2).OFFsetting = 122;
c1.capcontrol(2).ONsetting = 118;
c1.capcontrol(2).PTPhase = 3;
c1.capcontrol(2).PTRatio = 57.7;% was 60
c1.capcontrol(2).Vmax = 126;
c1.capcontrol(2).Vmin = 117;
c1.capcontrol(2).VoltOverride = 'True';
c1.capcontrol(2).Capacitor = c1.capacitor(1).Name;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c1.capcontrol(1) = [];% remove controls for the inactive capacitor
c1.capcontrol(2) = [];% remove controls for '46857CWT' since it always stays on.

%% change circuit settings
c1.circuit.basekv = '';
c1.circuit.Sequence = 'pos';
c1.circuit.pu = 1.01;%was 1.03
c1.circuit.bus1 = '';
c1.basevoltages = [115 12 0.208];%add 115kv as base

%% add "sourcebus" node to buslist
[AA,BB] = ismember('05465',c1.buslist.id);
coords = c1.buslist.coord(BB,:);
c1.buslist.id = ['SourceBus'; c1.buslist.id];
c1.buslist.coord = [coords; c1.buslist.coord];

%% add substation transformer
% if ~isfield(c1, 'transformer')
%    c1.transformer = dsstransformer('Name', 'Trafo'); 
% end
%     
c1.transformer(end+1) = c1.transformer(end);
c1.transformer(end).Name = 'TIMOTEO_1';
c1.transformer(end).Buses = {'SourceBus' '05465.1.2.3.0'};
c1.transformer(end).Conns = {'delta' 'wye'};
c1.transformer(end).kVs = [115 12];
c1.transformer(end).kVAs = [4200 4200];% changed from [10000 10000]
c1.transformer(end).XHL = 1.0871;
c1.transformer(end).sub = 'y';
c1.transformer(end).Rs = [0.103 0.103];

%% add regcontroller for substation transformer
n = dssregcontrol('Name',c1.transformer(end).Name);
n.transformer = c1.transformer(end).Name;
n.ptratio = 57.75;
n.EventLog = 'yes';
n.PTPhase = 3;
n.vlimit = 126;
n.revNeutral = [];
n.vreg = 122.8;% was 110
n.band = 3;
n.winding = 2;
c1.regcontrol = n;

%% Change load's kw and kvar according to Sunil's update
obj = excel2obj('C:\Work\Projects\2013\1665-CEC_Forecasting\Simulations\OpenDSS\Durox-Jens/Durox_SpotLoads_130730_Dec2012.xlsx');
obj = obj.SS_Custom_Spot_Load_Listings;
SpotNum = {obj(:,1).Spot_Number_}';% loads identifier from Sunil's document
kW = {obj(:,1).Total_kW_}';
kVar = {obj(:,1).Total_kVAR_}';
SpotNum = regexprep(SpotNum, '\$','_');SpotNum = regexprep(SpotNum, '\-','_');% Bring all names to uniform formatting
LoadNum = c1.load{:,1}.Name;% loads identifier from OpenDSS model
LoadNum = regexprep(LoadNum(:,1),'(\_\w+)','');% use only the first part of the string
LoadNum = regexprep(LoadNum,'l','');% remove the 'l' at the beginning of the string
SpotNum = regexprep(SpotNum(:,1),'(\_\w+)','');% use only the first part of the string
[aa,bb]=ismember(SpotNum, LoadNum);
[cc,dd]=ismember((1:1:129),bb);

for i=1:length(bb)% set the kw and kvar properties according to Sunil's document
    c1.load(bb(i)).Kw=cell2mat(kW(i));
    c1.load(bb(i)).Kvar=cell2mat(kVar(i));
end
c1.load(~cc)=[];% remove repeating loads from the load list

end