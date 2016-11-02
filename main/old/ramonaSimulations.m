clear;
close all;
clc;

%% Load circuit
c=load('data/f971_wPV.mat');c=c.c;
Days = struct('day',{'20121214', '20121218','20121219' },'weather',{'cloudy','overcast','clear'},...
    'date',{'5/25/2013','4/26/2013','6/2/2013'});

%% Add random PV
%_____________________________________________________________________________

%% get existing pv systems (exclude 500MW site)
pv0 = c.pvsystem;

pv_id = find([c.pvsystem.kVA] < 500);
pv = c.pvsystem(pv_id);
%% duplicate all existing system by 5 more times
pv8 = repmat(pv,1,8);
% rename system
for i = 1:length(pv8)
    pv8(i).Name = ['newPV' sprintf('_%03.0f',i)];
end

%% pick random bus (not including buses that already has pvsystems on them)
transbuslist ={};
for trans = 1:(length(c.transformer))
    transbuslist = {transbuslist{:}, c.transformer(trans).Buses{1},c.transformer(trans).Buses{2}};
end
adf={c.pvsystem.bus1}';
adf={adf{:},transbuslist{:}};
bus = setdiff(unique(cleanBus({c.line.bus1})'),adf); 
%% save new bus (DO NOT RUN THIS BLOCK since we don't want to regenerate random buses)
% bus_id = randperm(length(bus));
% bus_id = bus_id(1:length(pv8));
% save('data/f971busPV.mat','bus','bus_id');

%% replace PV buses with ramdom buses
load([normalizePath('$KLEISSLLAB24-1') 'database/gridIntegration/971/f971busPV.mat'])
pbcounter=0;
rdmpv=0;

for i = 1:length(pv8)
    if (~ismember(bus{bus_id(i)},cleanBus(transbuslist))) 
        pv8(i).bus1 = bus{bus_id(i)};
    else
        pv8(i).bus1 = '0909100C';
        pbcounter=pbcounter+1;
    end
end
disp([num2str(pbcounter) ' pb with transformers'])

%% add pv to original system
c.pvsystem = [pv0 pv8];
pv = [pv0 pv8];

save([pwd '/data/f971scenario2_pv.mat'],'pv');
%% reg settings
c.transformer(1).Name = ['t_' lower(c.transformer(1).Name)];
c.transformer(end).Name = lower(c.transformer(2).Name);
c.transformer(end).sub = 'y';
c.regcontrol = dssregcontrol;
c.regcontrol(1).transformer = lower(c.transformer(1).Name);
c.regcontrol(1).Name = c.regcontrol(1).transformer ;
c.regcontrol(1).winding = 2;
c.regcontrol(1).ptratio= 12000/120/sqrt(3);
c.regcontrol(1).vreg = 120.5;
c.regcontrol(1).band = 1.5;
c.regcontrol(2).transformer = c.transformer(end).Name;
c.regcontrol(2).Name = c.regcontrol(2).transformer ;
c.regcontrol(2).winding = 2;
c.regcontrol(2).ptratio= 12000/120/sqrt(3);
c.regcontrol(2).vreg = 121;
c.regcontrol(2).band = 1;
c.capacitor.Name = ['c' c.capacitor.Name];
c.capacitor.Numsteps=20;
c.capcontrol.Capacitor = c.capacitor.Name;
c.capcontrol.ONsetting=119;
c.capcontrol.OFFsetting=120;
c=rmfield(c,'capcontrol');
for l=1:length(c.capacitor)
    c.capacitor(l).Numsteps= 1;
end
% Energymeter.
c.energymeter = dssenergymeter;
c.energymeter.Name='substation';
c.energymeter.element=['transformer.' c.transformer(end).Name];
save('data/f971_test.mat','c')

%% add PV profils
for iii=1:3
% iii=3
    circuitPath = [pwd '/data/f971_test.mat']; % file generated upper in the page
    pvForecastPath = [normalizePath('$KLEISSLLAB24-1') 'database/gridIntegration/971/' Days(iii).day]; 
    desag_ag = [1 1]; 
    if (desag_ag==0); profil_type = 'aggregated'; else profil_type = 'disaggregated'; end
    outPath = ['tmp/'];
    addBigPlant=[0 0];
    plotOption = 0;

    fprintf(['Creating ' profil_type ' pv profil:\n'])
    [c p] = AddPVForecast(circuitPath, pvForecastPath, desag_ag, outPath, addBigPlant, plotOption);
    %% change to generator
    if isfield(c,'generator'),offset=length(c.generator);else offset=0,c.generator(1) = dssgenerator; end
    for i = 1:length(c.pvsystem)
        % Panel kW = Pmpp (in kW @1kW/m2 and 25 C) * Irradiance (in kW/m2) * Factor(@actual T)
        c.generator(i+offset).Name = ['gen_' num2str(i)];
        c.generator(i+offset).Kw = c.pvsystem(i).kVA;
        c.generator(i+offset).Kvar = 0;
        c.generator(i+offset).kv = 12;
        c.generator(i+offset).phases = 1;
        c.generator(i+offset).Vmaxpu = 1.5;
        c.generator(i+offset).Vminpu = 0.6;
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
    %% add Loadshape
%     date={'5/25/2013','4/26/2013','6/2/2013'}
    [ c  p] = AddLS( c , pwd);%, date{iii} );
    %  c.loadshape(2).Mult=0:1/2880:1-1/2880;

    %% TESTING
    summprod=0;
    for i =2:length(c.generator)
        summprod = summprod + c.generator(i).Kw*c.loadshape(find(ismember({c.loadshape.Name},c.generator(i).Daily))).Mult(1500);
        prof1500(i-1) = c.loadshape(find(ismember({c.loadshape.Name},c.generator(i).Daily))).Mult(1500);
        yop(i-1) = isequal(c.loadshape(find(ismember({c.loadshape.Name},c.generator(i).Daily))).Mult, zeros(1,2880)+0.0001);
    end
    summprod
    c_bk=c;

    %% Simulation no pv
%     c= load('data/f520_test');c=c.c;
    cbk=c;
    for k=1:length(c.generator)
        c.generator(k).Enabled='false';
    end
    disp('No PV:')
    c.dataNoPV =dssSimulation(c,'daily',30);
    
    %% Simulation 4% (Real)
    for k=1:43
        c.generator(k).Enabled='true';
    end
    for k=44:length(c.generator)
        c.generator(k).Enabled='false';
    end
%     sum([c.generator{1:43}.Kw])/sum([c.load.Kw])
    disp('4% PV pen.:')
    c.dataPV4 =dssSimulation(c,'daily',30);
    %% Simulation 18%
    for k=1:196
        c.generator(k).Enabled='true';
    end
    for k=197:length(c.generator)
        c.generator(k).Enabled='false';
    end
%     sum([c.generator{1:100}.Kw])/sum([c.load.Kw])
    disp('18% PV pen.:')
    c.dataPV18 =dssSimulation(c,'daily',30);
    
    %% Simulation 30%
    for k=1:324
        c.generator(k).Enabled='true';
    end
    for k=325:length(c.generator)
        c.generator(k).Enabled='false';
    end
%     sum([c.generator{1:324}.Kw])/sum([c.load.Kw])
    disp('30% PV pen.:')
    c.dataPV30 =dssSimulation(c,'daily',30);
    
    %% Simulation 50%
    c.generator=cbk.generator;
    sFactor=1.39;
    for k=1:length(c.generator)
        c.generator(k).Enabled='true';
    end
    for k=1:length(c.generator)
        c.generator(k).Kw=c.generator(k).Kw*sFactor;
    end
%     sum([c.generator.Kw])/sum([c.load.Kw])
    disp('50% PV pen.:')
    c.dataPV50 =dssSimulation(c,'daily',30);

    %% Simulation 75%
    c.generator=cbk.generator;
    sFactor=2.075;
    for k=1:length(c.generator)
        c.generator(k).Enabled='true';
    end
    for k=1:length(c.generator)
        c.generator(k).Kw=c.generator(k).Kw*sFactor;
    end
%     sum([c.generator.Kw])/sum([c.load.Kw])
    disp('75% PV pen.:')
    c.dataPV75 =dssSimulation(c,'daily',30);
    
    %% Simulation 100%
    c.generator=cbk.generator;
    sFactor=2.77;
    for k=1:length(c.generator)
        c.generator(k).Enabled='true';
    end
    for k=1:length(c.generator)
        c.generator(k).Kw=c.generator(k).Kw*sFactor;
    end
%     sum([c.generator.Kw])/sum([c.load.Kw])
    disp('100% PV pen.:')
    c.dataPV100 =dssSimulation(c,'daily',30);
%     save('data/f520Solved','c');
    
    namelist={'dataPV4','dataPV18','dataPV30','dataPV50','dataPV75','dataPV100'};
    
    for l=1:6 % how many cases
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
        mkdir([pwd '/tmp/f971_' namelist{l}(7:end) '/'])
        result_dir=[pwd '/tmp/f971_' namelist{l}(7:end) '/'];
        DailyPlot(c.(namelist{l}), c.dataNoPV, result_dir,...
            ['f971_' namelist{l}(7:end) '_' Days(iii).weather],1,[13 13.5 14])
    end
    save([pwd '/data/f971_SimResults_' Days(iii).weather '.mat'],'c');
end
%% Plot comparison between all cases
clear Data;
Data ={};
PVpen= [0, 4, 18, 30, 50, 75, 100];
for w=1:3   %weather 
    resu_file=[pwd '/data/f971_SimResults_' Days(w).weather '.mat'];
    c=load(resu_file);c=c.c;
    Data{w,1} = c.dataNoPV;
    for g=2:length(PVpen) 
        Data{w,g} = c.(['dataPV' num2str(PVpen(g))]);
    end
end


leg={'cloudy','overcast','clear'};
mkdir([pwd 'tmp/f971'])
outPath = 'tmp/f971';
dataplotPVpen( PVpen, Data ,leg, outPath)