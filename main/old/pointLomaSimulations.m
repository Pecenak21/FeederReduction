
%% Summary:
% 1. output: data/f480sc1c1.mat -> circuit with original PV systems (85)
% 2. output: data/f480sc1c2.mat -> circuit with 340 pv systems 
% 3. run 5 cases: 1->5 includes diff PV penetration (ori 15%, 30, 60, 43, 50) and 6 is with big pv system (ori 1 at 16 different locations  

%%
clear;
%% Scenario 1:

Days = struct('day',{'20121214', '20121218','20121219' },'weather',{'cloudy','overcast','clear'},...
    'date',{'12/14/2012','12/18/2012','12/19/2012'});

% Load a circuit with pv systems with azi and tilt angles
load('data/f480_wPV.mat');

load data/CABRILLO_pvs_Ti_Az;
if length(pvs)==85;
    pvs=[pvs , pvs, pvs, pvs];
    save([pwd '/data/' c.circuit.Name '_pvs_Ti_Az.mat'], 'pvs');
end
% %% PV sytems
for bla=1:length(c.pvsystem)
    c.pvsystem(bla).kvar = 0;
end

% %% Update sourcebus
c.buslist.coord(end+1,:) = c.buslist.coord(3,:);c.buslist.id(end+1) = {'SourceBus'};

% %% modify transformer
c.transformer.kVs = [12 12];

% %% Add new transformer and regcontrol objects
% needed to know how phases are in each transformer
lineE = regexp(c.line(:).bus1,'\.','split');
aa = lineE(:);
for i=1:length(lineE)
    aaa(i) = aa{i}(1);
end
% c.basevoltages = [12 12 12];
c.circuit.bus1 = ''; 
c.transformer(2) = dsstransformer;
c.transformer(2).Name = 'cabrillo';
c.transformer(2).Buses = {'SourceBus' '0480'};% '04801.1.2.3.0'};%{'048030' '048031'};
c.transformer(2).Phases = 3;
c.transformer(2).kVAs = [10000 10000];
c.transformer(2).kVs = [12 12];
c.transformer(2).Conns = {'wye' 'wye'};
c.transformer(2).sub = 'yes';
c.transformer(2).Rs = [0.103       0.103];
c.transformer(2).Windings = 2;
c.transformer(2).XHL = 1.0871;
c.transformer(2).Wdg = 2;
c.transformer(2).MaxTap = 1.1;
c.transformer(2).MinTap = 0.9;

c.transformer(3) = dsstransformer;
c.transformer(3).Name = 'T3';
c.transformer(3).Buses = {'048043' '048044'};%{'04804707' '04804711'};
c.transformer(3).Phases = c.line(find(ismember(aaa,c.transformer(3).Buses{1}))).Phases;
c.transformer(3).Windings = 2;
c.transformer(3).kVAs = [10000 10000];
c.regcontrol(1) = dssregcontrol;
c.regcontrol(1).Name = 't3';
c.regcontrol(1).transformer = 'T3';
c.regcontrol(1).PTPhase = c.transformer(3).Phases;
c.regcontrol(1).winding = 2;
c.regcontrol(1).ptratio = 12*1000/120/sqrt(3);
c.regcontrol(1).delay = 15;

c.transformer(4) = dsstransformer;
c.transformer(4).Name = 'T4';
c.transformer(4).Buses = {'048056','048057'}; %{'048045' '048046'};%{'048082' '048083'};
c.transformer(4).Phases = c.line(find(ismember(aaa,c.transformer(4).Buses{1}))).Phases;
c.transformer(4).Windings = 2;
c.transformer(4).kVAs = [8000 8000];
c.regcontrol(2) = dssregcontrol;
c.regcontrol(2).Name = 't4';
c.regcontrol(2).transformer = 'T4';
c.regcontrol(2).PTPhase = c.transformer(4).Phases;
c.regcontrol(2).winding = 2;
c.regcontrol(2).ptratio = 12*1000/120/sqrt(3);
c.regcontrol(2).delay = 15+10;


for i=3:length(c.transformer)
%     if (c.transformer(i).Phases > 1)
        c.transformer(i).kVs = [12 12];
%     else
%         c.transformer(i).kVs = [12/sqrt(3) 12/sqrt(3)];
%     end
%     c.transformer(i).Conns = {'wye', 'wye'};
    c.transformer(i).Wdg = 2; % to change the secondary winding
    c.transformer(i).MaxTap = 1.1;
    c.transformer(i).MinTap = 0.9;
    c.regcontrol(i-2).band = 1.5;%7.2
    c.regcontrol(i-2).vreg = 120.5;
    c.regcontrol(i-2).ptratio = 12*1000/120/sqrt(3);
    c.regcontrol(i-2).vlimit = 126;
end
c.regcontrol(1).vreg = 120.5;
% c.transformer(3).kVs = [12 12*(1+0.02)];
c.regcontrol(3) = dssregcontrol;
c.regcontrol(3).Name = c.transformer(2).Name;
c.regcontrol(3).transformer = c.transformer(2).Name;
c.regcontrol(3).PTPhase = c.transformer(3).Phases;
c.regcontrol(3).winding = 2;
c.regcontrol(3).ptratio = 12*1000/120/sqrt(3);
% c.regcontrol(3).delay = 10;
c.regcontrol(3).band = 1;%7.2
c.regcontrol(3).vreg = 121;

% %% Add capcontrol objects at each variable capacitor bank
LINE = c.line;

for j=1:length(c.line)
    bb=regexp([c.line(j).bus1], '\.', 'split');
    LINE(j).bus1=bb(1);
    bb=regexp([c.line(j).bus2], '\.', 'split');
    LINE(j).bus2=bb(1);
end
% Create capcontrol
j=[3,4,5,6];
c.capcontrol = dsscapcontrol;
for cpctr=1:4
    c.capcontrol(cpctr) = dsscapcontrol;
    c.capcontrol(cpctr).Name = lower(c.capacitor(j(cpctr)).Name);
    AA = ismember([LINE(:).bus1]',c.capacitor(j(cpctr)).Bus1);
    if (sum(AA)==0)
        AA = ismember([LINE(:).bus2]',c.capacitor(j(cpctr)).Bus1);
    end
    if (length(c.line(AA).Name)>1) && (length(c.line(AA).Name)<4)
        c.capcontrol(cpctr).Element = ['line.' c.line(AA).Name{1}];
    end
    if (length(c.line(AA).Name)>3)
        c.capcontrol(cpctr).Element = ['line.' c.line(AA).Name];
    end
    c.capcontrol(cpctr).Capacitor = lower(c.capacitor(j(cpctr)).Name);
end
% change settings
for j=1:length(c.capcontrol)
    c.capcontrol(j).Type = 'voltage';
    c.capcontrol(j).OFFsetting = 120;
    c.capcontrol(j).ONsetting = 118;
    c.capcontrol(j).PTRatio = 12*1000/120/sqrt(3);
    c.capcontrol(j).VoltOverride = 'yes';
    c.capcontrol(j).Vmax = 126;
    c.capcontrol(j).Vmin = 114;
    c.capcontrol(j).CTRatio=300;
end
% %% Correct pv phase and bus
LINE = c.line;

for j=1:length(c.line)
    bb=regexp([c.line(j).bus1], '\.', 'split');
    LINE(j).bus1=bb(1);
    bb=regexp([c.line(j).bus2], '\.', 'split');
    LINE(j).bus2=bb(1);
end
for kk=1:length(c.pvsystem)
    [ab abc] = ismember(c.pvsystem(kk).bus1,{LINE.bus1});
    if abc==0;
        [ab abc] = ismember(c.pvsystem(kk).bus1,{LINE.bus2});
        c.pvsystem(kk).bus1 = c.line(abc).bus2;
    else
        c.pvsystem(kk).bus1 = c.line(abc).bus1;
    end
    c.pvsystem(kk).phases = c.line(abc).Phases;
    if (c.line(abc).Phases>3) | (c.line(abc).Phases<1)
        c.line(abc).Phases
    end
end
% %% Save 
run('def_addpath.m')
save([pwd '/data/f480sc1c1.mat'],'c','p','glc','d'); % initial PV System

% %% Add random PV
%_____________________________________________________________________________

% %% get existing pv systems (exclude 2MW site)
pv0 = c.pvsystem;

pv_id = find([c.pvsystem.kVA] < 990);
pv = c.pvsystem(pv_id);
% %% duplicate all existing system by 4  more times
pv3 = repmat(pv,1,3);
% rename system
for i = 1:length(pv3)
    pv3(i).Name = ['newPV' sprintf('_%03.0f',i)];
end

% %% pick random bus (not including buses that already has pvsystems on them)
bus = setdiff(unique({LINE.bus1}'),{c.pvsystem.bus1}'); 
% %% save new bus (DO NOT RUN THIS BLOCK since we don't want to regenerate random buses)
% bus_id = randperm(length(bus));
% bus_id = bus_id(1:length(pv3));
% save('data/f480busPV_1.mat','bus','bus_id');
% %% replace PV buses with ramdom buses
load data/f480busPV.mat
cba = [122, 135, 142, 143, 139, 141, 130, 131, 125, 126,8502,8503,8508,8516,8517,8207,8208,8234,8235,8203,8202];
pbcounter=0;
rdmpv=0;
transbuslist ={};
for trans = 1:(length(c.transformer)-1)
    transbuslist = {transbuslist{:}, c.transformer(trans+1).Buses{1},c.transformer(trans+1).Buses{2}};
end
for i = 1:length(pv3)
    if (~ismember(bus{bus_id(i)},transbuslist)) 
%             & (str2num(bus{bus_id(i)}(1:5))>=4804) & (~eq(str2num(bus{bus_id(i)}(1:5)),4806))...
%             & (~eq(str2num(bus{bus_id(i)}(1:5)),4807))& (~eq(str2num(bus{bus_id(i)}(1:5)),4808))
        pv3(i).bus1 = bus{bus_id(i)};
    else
        pbcounter=pbcounter+1;
        rdmpv=rdmpv+1;
        if rdmpv>length(cba)
            rdmpv=rdmpv-length(cba);
        end
        disp(['bus ' bus{bus_id(i)} ' pv ' pv3(i).Name '\n' num2str(pv3(i).kVA) ])
        pv3(i).bus1 = ['0480' num2str(cba(rdmpv))];
    end
end
disp([num2str(pbcounter) ' pb with transformers'])

% %% add pv to original system
c.pvsystem = [pv0 pv3];
pv = [pv0 pv3];

%save([pwd '/data/f480scenario2_pv.mat'],'pv');

% %% fix the phases for PVsystems
pv = load('data/f480scenario2_pv.mat'); pv = pv.pv; 
pv2 = pv(1:170);

% %% save pv system
def_addpath;
save([pwd '/data/scenario2_pv.mat'],'pv');
c.pvsystem = pv;
% correct pv phase
LINE = c.line;

for j=1:length(c.line)
    bb=regexp([c.line(j).bus1], '\.', 'split');
    LINE(j).bus1=bb(1);
    bb=regexp([c.line(j).bus2], '\.', 'split');
    LINE(j).bus2=bb(1);
end
for kk=1:length(c.pvsystem)
    pvbus = regexp(c.pvsystem(kk).bus1,'\.','split');
    c.pvsystem(kk).bus1 = pvbus(1);
end
for kk=1:length(c.pvsystem)
    [ab abc] = ismember(c.pvsystem(kk).bus1,{LINE.bus1});
    if abc==0;
        [ab abc] = ismember(c.pvsystem(kk).bus1,{LINE.bus2});
        c.pvsystem(kk).bus1 = c.line(abc).bus2;
    else
        c.pvsystem(kk).bus1 = c.line(abc).bus1;
    end
    c.pvsystem(kk).phases = c.line(abc).Phases;
    if (c.line(abc).Phases>3) | (c.line(abc).Phases<1)
        c.line(abc).Phases
    end
end
c=rmfield(c,'capcontrol');
for l=1:length(c.capacitor)
    c.capacitor(l).Numsteps= 1;
end
% Energymeter.
c.energymeter = dssenergymeter;
c.energymeter.Name='substation';
c.energymeter.element=['transformer.' c.circuit.Name];
save([pwd '/data/f480sim.mat'],'c','glc');

for iii=1:3
    
    circuitPath = ['data/f480sim.mat']; % file generated upper in the page
    pvForecastPath = [normalizePath('$KLEISSLLAB24-1') 'database/gridIntegration/480/USI_1_2/' Days(iii).day]; ;
    desag_ag = [1 1]; 
    if (desag_ag==0); profil_type = 'aggregated'; else profil_type = 'disaggregated'; end
    outPath = ['tmp/f480'];
    addBigPlant=[0 0];
    plotOption = 0;

% %% Add pv profil

    fprintf(['Creating ' profil_type ' pv profil:\n'])
     [c p] = AddPVForecast480(circuitPath, pvForecastPath, desag_ag, outPath,addBigPlant, plotOption);
    for i=1:length(c.loadshape)
        if max([c.loadshape(i).Mult])<0.001;
            test(i)=0;
        else test(i)=1;
        end
    end

    offset=0;
    for i = 1:length(c.pvsystem)
        % Panel kW = Pmpp (in kW @1kW/m2 and 25 C) * Irradiance (in kW/m2) * Factor(@actual T)
        c.generator(i+offset)=dssgenerator;
        c.generator(i+offset).Name = ['gen_' num2str(i)];
        c.generator(i+offset).Kw = c.pvsystem(i).kVA;
        c.generator(i+offset).Kvar = 0;
        c.generator(i+offset).kv = 12;
        c.generator(i+offset).phases = 1;
        c.generator(i+offset).Vmaxpu = 1.8;
        c.generator(i+offset).Vminpu = 0.3;
        c.generator(i+offset).Pf=1;
%         c.generator(i+offset).Debugtrace='yes';
        if c.generator(i).Kw>60
            c.generator(i+offset).bus1 = [cleanBus(c.pvsystem(i).bus1) '.1.2.3'];
        else 
            c.generator(i+offset).bus1 = c.pvsystem(i).bus1;
        end
        c.generator(i+offset).Daily = c.pvsystem(i).daily;
    end
    c=rmfield(c,'pvsystem');
    for k=1:length(c.generator)
        c.generator(k).Enabled='true';
    end
    
%     % Add loadshape
    date=Days(iii).date;
    [ c  p] = AddLS( c , outPath, date );
    
        c.transformer(1).Name = ['t' c.transformer(1).Name];
    for x=1:length(c.transformer)
        c.transformer(x).Name = lower(c.transformer(x).Name);
    end
    for x=1:length(c.regcontrol)
        c.regcontrol(x).Name = lower(c.regcontrol(x).Name);
        c.regcontrol(x).transformer = lower(c.regcontrol(x).transformer);
    end
    cbk=c;
    %%
    cbk=c;
    for k=1:length(c.generator)
        c.generator(k).Enabled='false';
    end
    disp('No PV:')
    c.dataNoPV =dssSimulation(c,'daily',30);
    
    for k=1:85
        c.generator(k).Enabled='true';
    end
    for k=86:length(c.generator)
        c.generator(k).Enabled='false';
    end
    sum([c.generator{1:85}.Kw])/sum([c.load.Kw])
    disp('16% PV pen.:')
    c.dataPV16 =dssSimulation(c,'daily',30);
    %
    for k=1:length(c.generator)
        c.generator(k).Enabled='true';
    end
    sum([c.generator.Kw])/sum([c.load.Kw])%61
    disp('61% PV pen.:')
    c.dataPV61 =dssSimulation(c,'daily',30);
    
    for k=1:length(c.generator)
        c.generator(k).Enabled='true';
    end
    for k=256:340
        c.generator(k).Enabled='false';
    end
    sum([c.generator{1:255}.Kw])/sum([c.load.Kw])% 46
    disp('46% PV pen.:')
    c.dataPV46 =dssSimulation(c,'daily',30);
    
%     c.generator=cbk.generator;
    x=0;
    for k=1:85
        c.generator(k).Enabled='true';
    end
    for k=85:340
        if x==1
%             summ = summ + c.pvsystem(k).kVA;
            c.generator(k).Enabled='true';
            x=x-1;
        else
            x=x+1;
            c.generator(k).Enabled='false';
        end
    end
    sum([c.generator{ismember({c.generator.Enabled},'true')}.Kw])/sum([c.load.Kw])%37
    disp('37% PV pen.:')
    c.dataPV37 =dssSimulation(c,'daily',30);
    
    x=0;
    for k=86:251
        if x==1
            c.generator(k).Enabled='true';
            x=x-1;
%             counter=counter+1;
        else
            c.generator(k).Enabled='false';
            x=x+1;
        end
    end
    for k=[252:260,262:340]
        c.generator(k).Enabled='false';
    end
    sum([c.generator{ismember({c.generator.Enabled},'true')}.Kw])/sum([c.load.Kw])%30
    disp('30% PV pen.:')
    c.dataPV30 =dssSimulation(c,'daily',30);
    

    c.generator=cbk.generator;
    sFactor=1.24;
    for k=1:length(c.generator)
        c.generator(k).Enabled='true';
        c.generator(k).Kw=c.generator(k).Kw*sFactor;
    end
    sum([c.generator{ismember({c.generator.Enabled},'true')}.Kw])/sum([c.load.Kw])
    disp('75% PV pen.:')
    c.dataPV75 =dssSimulation(c,'daily',30);
    
    
    c.generator=cbk.generator;
    sFactor=1.64;
    for k=1:length(c.generator)
        c.generator(k).Enabled='true';
        c.generator(k).Kw=c.generator(k).Kw*sFactor;
    end
    sum([c.generator{ismember({c.generator.Enabled},'true')}.Kw])/sum([c.load.Kw])
%     c=addMonitors( c, {'generator'});
    disp('100% PV pen.:')
    c.dataPV100 =dssSimulation(c,'daily',30);
    summ=0;
    for o=1:length(c.generator)
        summ=summ+c.generator(o).Kw*c.loadshape(find(ismember({c.loadshape.Name},c.generator(o).Daily))).Mult(1529);
    end
    summ;
    namelist={'dataPV16','dataPV30','dataPV37','dataPV46','dataPV61','dataPV75','dataPV100'};
    
    for l=1:7
        sunset = 17; % 5pm 
%         c2pnames= fieldnames(c.(namelist{l}).Cap2Plot);
%         for k=1:length(c.(namelist{l}).Cap2Plot)
%             c.(namelist{l}).Cap2Plot.(c2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Cap2Plot.(c2pnames{k})(2880*sunset/24:end);
%         end
        t2pnames= fieldnames(c.(namelist{l}).Tap2Plot);
        for k=1:length(c.(namelist{l}).Tap2Plot)
            c.(namelist{l}).Tap2Plot.(t2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Tap2Plot.(t2pnames{k})(2880*sunset/24:end);
        end
        c.(namelist{l}).LineLoss(2880*sunset/24:end,:)=c.dataNoPV.LineLoss(2880*sunset/24:end,:);
        c.(namelist{l}).TotalLoss(2880*sunset/24:end,:)=c.dataNoPV.TotalLoss(2880*sunset/24:end,:);
        c.(namelist{l}).TotalPower(2880*sunset/24:end,:)=c.dataNoPV.TotalPower(2880*sunset/24:end,:);
        c.(namelist{l}).VoltMaxMin(2880*sunset/24:end,:)=c.dataNoPV.VoltMaxMin(2880*sunset/24:end,:);
        c.(namelist{l}).Voltage(2880*sunset/24:end,ismember(c.(namelist{l}).nodeName,c.dataNoPV.nodeName)) = ...
            c.(namelist{l}).Voltage(2880*sunset/24:end,ismember(c.(namelist{l}).nodeName,c.dataNoPV.nodeName));
        c.(namelist{l}).Capcontrol(2880*sunset/24:end,:) = c.dataNoPV.Capcontrol(2880*sunset/24:end,:);
        c.(namelist{l}).Regulation(2880*sunset/24:end,:) = c.dataNoPV.Regulation(2880*sunset/24:end,:);
        mkdir([pwd '/tmp/f480pen' namelist{l}(7:end) '/'])
        result_dir=[pwd '/tmp/f480pen' namelist{l}(7:end) '/'];
        DailyPlot(c.(namelist{l}), c.dataNoPV, result_dir,...
            ['f480' namelist{l}(7:end) '_' Days(iii).weather],1,[13 13.5 14])
    end
    save([pwd '/data/f480_SimResults_' Days(iii).weather '.mat'],'c');
end

%% Plot comparison between all cases
clear Data;

PVpen= [0, 16, 30, 37, 46, 61, 75, 100];
Data ={};
for w=[1,3,2]  %weather 
    c=load([pwd '/data/f480_SimResults_' Days(w).weather '.mat']);c=c.c;
    for g=1:7 % number of cases
        if g==1,
            Data{w,g} = c.dataNoPV;
        end
        Data{w,g+1} = c.(['dataPV' num2str(PVpen(g+1))]);
    end
end


leg={'cloudy','clear','overcast'};
outPath = 'tmp/f480';
dataplotPVpen( PVpen, Data ,leg, outPath)

