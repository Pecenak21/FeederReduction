function [c1] = circuit_changes_Centaur(c1, node_bk)

LineDSS = 'SCE_Centaur_missing_lines.dss';

%% line impedances need to have tiny values rather than 0s
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
n.element = 'Line.03265_GS1305_1_03265_1_03265';
n.terminal = 1;
c1.energymeter = n;

%% add missing transformers
for i=1:length(c1.transformer)
    c1.transformer(i).kVs = [12 0.208];
    c1.transformer(i).Windings = 2;
end
n = c1.transformer(end);
n.Name = 't_707';
n.Buses = {'664.1.2.3', '690.1.2.3'};
c1.transformer(end+1) = n;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 't_708';
n.Buses = {'664.1.2.3', '699.1.2.3'};
c1.transformer(end+1) = n;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 't_712';
n.Buses = {'664.1.2.3', '693.1.2.3'};
c1.transformer(end+1) = n;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 't_729';
n.Buses = {'664.1.2.3', '726.1.2.3'};
c1.transformer(end+1) = n;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 't_730';
n.Buses = {'664.1.2.3', '727.1.2.3'};
c1.transformer(end+1) = n;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 't_755';
n.Buses = {'745.1.2.3', '741.1.2.3'};
c1.transformer(end+1) = n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 't_758';
n.Buses = {'751.1.2.3', '747.1.2.3'};
c1.transformer(end+1) = n;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 't_759';
n.Buses = {'737.1.2.3', '734.1.2.3'};
c1.transformer(end+1) = n;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 't_765';
n.Buses = {'737.1.2.3', '738.1.2.3'};
c1.transformer(end+1) = n;
% 
% % add transformer bus coordinates
for i=1:length(c1.transformer)
    bus1 = c1.transformer(i).Buses(1,1);
    bus = regexprep(bus1,'(\.\d+)+$','');
    if ~ismember(bus,c1.buslist.id)
        [AA,BB]=ismember(bus,node_bk.NodeID);
        if BB > 0
            c1.buslist.id = [bus;c1.buslist.id];
            coordX = str2double(node_bk.CoordX(BB,:));
            coordY = str2double(node_bk.CoordY(BB,:));
            coords = [coordX,coordY];
            c1.buslist.coord = [coords;c1.buslist.coord];
        end
    end  
    bus2 = c1.transformer(i).Buses(1,2);
    bus = regexprep(bus2,'(\.\d+)+$','');
    if ~ismember(bus,c1.buslist.id)
        [AA,BB]=ismember(bus,node_bk.NodeID);
        if BB > 0
            c1.buslist.id = [bus;c1.buslist.id];
            coordX = str2double(node_bk.CoordX(BB,:));
            coordY = str2double(node_bk.CoordY(BB,:));
            coords = [coordX,coordY];
            c1.buslist.coord = [coords;c1.buslist.coord];
        end
    end 
end

%% add regcontroller for pv transformers
n = dssregcontrol('Name',c1.transformer(1).Name);
n.transformer = c1.transformer(1).Name;
n.vreg = 122.8;
n.band = 3;
n.ptratio = 57.7;
n.EventLog = 'yes';
n.PTPhase = 3;
n.vlimit = 126;
n.revNeutral = [];
n.winding = 2;
c1.regcontrol = n;
for i=2:length(c1.transformer)
    n.Name = c1.transformer(i).Name;
    n.transformer = c1.transformer(i).Name;
    c1.regcontrol(end+1) = n;
end

%% add remaining PV manually
PV = dsspvsystem('Name', '802');
PV.bus1 = '784';
PV.phases = 3;
PV.kv = 0.208;
PV.kVA = 399;
PV.pf = 1;
PV.EffCurve = 'eff';
PV.PTCurve = 'PvsT';
PV.Tdaily = 'Temp';
PV.Pmpp = 399;
PV.cutin = 0.01;
PV.cutout = 1.2;
PV.conn = 'delta';
PV.irradiance = 1;
PV.temperature = 25;
c1.pvsystem = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '807';
PV.bus1 = '780';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '809';
PV.bus1 = '793';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '815';
PV.bus1 = '787';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '709';
PV.bus1 = '690';
PV.kVA = 238;
PV.Pmpp = 238;
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '703';
PV.bus1 = '699';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '705';
PV.bus1 = '693';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '731';
PV.bus1 = '726';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '728';
PV.bus1 = '727';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '769';
PV.bus1 = '741';
PV.kVA = 191;
PV.Pmpp = 191;
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '763';
PV.bus1 = '747';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '761';
PV.bus1 = '734';
c1.pvsystem(end+1) = PV;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PV.Name = '756';
PV.bus1 = '738';
c1.pvsystem(end+1) = PV;
% 
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
temp.interval = 1;
temp.temp = [25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25	25];
c1.tshape = temp;

%% remove additional lines (during conversion 4 switches have been added that need to be replaced with PVs)
PVs=c1.pvsystem(:).Name;
[aa,bb]=ismember(PVs, c1.line(:).Name);
bb=bb(aa);
c1.line(bb)='';

%% Add missing lines
line = dssparse(LineDSS);%missing lines are stored in a dss file. that appeared to be the simplest way to get them into the circuit
for i=1:length(line.line)
    n = line.line(i);
    c1.line(end+1)=n;
    bus1 = regexprep(n.bus1,'(\.\d+)+$','');
    bus2 = regexprep(n.bus2,'(\.\d+)+$','');
    if ~ismember(bus1,c1.buslist.id)% add bus coordinates
        [AA,BB]=ismember(bus1,node_bk.NodeID);
        if BB > 0
            c1.buslist.id = [bus1;c1.buslist.id];
            coordX = str2double(node_bk.CoordX(BB,:));
            coordY = str2double(node_bk.CoordY(BB,:));
            coords = [coordX,coordY];
            c1.buslist.coord = [coords;c1.buslist.coord];
        end
    end
    if ~ismember(bus2,c1.buslist.id)
        [AA,BB]=ismember(bus2,node_bk.NodeID);
        if BB > 0
            c1.buslist.id = [bus2;c1.buslist.id];
            coordX = str2double(node_bk.CoordX(BB,:));
            coordY = str2double(node_bk.CoordY(BB,:));
            coords = [coordX,coordY];
            c1.buslist.coord = [coords;c1.buslist.coord];
        end
    end
end

%% remove additional lines (during conversion all loads were connected to the circuit through switches. These
% are not required.
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

%% add missing loads
% n = c1.load(end);
% n.Name = 'lP5611340_T_03265';
% n.Phases = 3;
% n.bus1 = 'P5611340_03265';
% n.Kw = 0;
% n.Kvar = 0;
% c1.load(end+1) = n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = c1.load(end);
n.Name = 'lP5541088_T_03265';
n.Phases = 3;
n.bus1 = 'P5541088_03265';
n.Kw = 93.39081;
n.Kvar = 40.90286446;
c1.load(end+1) = n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 'lP5565162_T_03265';
n.Phases = 3;
n.bus1 = 'P5565162_03265';
n.Kw = 280.183638;
n.Kvar = 122.7135022;
c1.load(end+1) = n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.Name = 'lP5601881_T_03265';
n.Phases = 3;
n.bus1 = 'P5601881_1_T_03265';
n.Kw = 56.036727;
n.Kvar = 24.54270018;
c1.load(end+1) = n;

%% add voltage setting to loads
for i=1:length(c1.load)
    if (c1.load(i).Phases > 1)
        c1.load(i).Kv = 12;
    else
        c1.load(i).Kv = 12/sqrt(3);
    end
    c1.load(i).NumCust = 1;
end

%% change loads model to see if that affects something
% for i=1:length(c1.load)
% %     c1.load(i).Model = 1;
%     c1.load(i).ZIPV = [0.7           0         0.3         0.7
%     0         0.3         0.5];
% end

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
c1.capacitor(2).Kv = 12*0.95;
c1.capacitor(2).Name = ['cap_' c1.capacitor(2).Name];
c1.capacitor(3).Kv = 12*0.95;
c1.capacitor(3).Name = ['cap_' c1.capacitor(3).Name];
c1.capacitor(1) = [];% remove inactive capacitor

%% change capcontrol settings
c1.capcontrol(2).Type = 'voltage';
c1.capcontrol(2).CTRatio = 100;
c1.capcontrol(2).Eventlog = 'True';
c1.capcontrol(2).OFFsetting = 122.5;
c1.capcontrol(2).ONsetting = 118;
c1.capcontrol(2).PTPhase = 3;
c1.capcontrol(2).PTRatio = 57.7;
c1.capcontrol(2).Vmax = 126;
c1.capcontrol(2).Vmin = 117;
c1.capcontrol(2).VoltOverride = 'True';
c1.capcontrol(2).Capacitor = c1.capacitor(1).Name;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c1.capcontrol(1) = [];% remove controls for the inactive capacitor
c1.capcontrol(2) = [];% remove controls for '2210560E' since it always stays on.

%% change circuit settings
c1.circuit.basekv = '';
c1.circuit.Sequence = 'pos';
c1.circuit.pu = 1.01;
c1.circuit.bus1 = '';
c1.basevoltages = [115 12 0.208];%add 115kv as base

%% add "sourcebus" node to buslist
[AA,BB] = ismember('03265',c1.buslist.id);
coords = c1.buslist.coord(BB,:);
c1.buslist.id = ['SourceBus'; c1.buslist.id];
c1.buslist.coord = [coords; c1.buslist.coord];

%% add substation transformer
c1.transformer(end+1) = c1.transformer(end);
c1.transformer(end).Name = 'SAN_BERDO';
c1.transformer(end).Buses = {'SourceBus' '03265.1.2.3.0'};
c1.transformer(end).Conns = {'delta' 'wye'};
c1.transformer(end).kVs = [115 12];
c1.transformer(end).kVAs = [5700 5700];% changed from [10000 10000]
c1.transformer(end).XHL = 1.0871;
c1.transformer(end).sub = 'y';
c1.transformer(end).Rs = [0.103 0.103];

%% add regcontroller for substation transformer
n = dssregcontrol('Name',c1.transformer(end).Name);
n.transformer = c1.transformer(end).Name;
n.ptratio = 57.7;
n.EventLog = 'yes';
n.PTPhase = 3;
n.vlimit = 126;
n.revNeutral = [];
n.vreg = 122.8;
n.band = 3;
n.winding = 2;
c1.regcontrol(end+1) = n;
% c1.regcontrol = n;

end