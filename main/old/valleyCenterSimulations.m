clear;
close all;
clc;

%% Load circuit
c=load('data/f909_wPV.mat');c=c.c;
Days = struct('day',{'20121214', '20121218','20121219' },'weather',{'cloudy','overcast','clear'},...
    'date',{'5/25/2013','4/26/2013','6/2/2013'});

%% Add random PV
%_____________________________________________________________________________

%% get existing pv systems (exclude 500MW site)
pv0 = c.pvsystem;

pv_id = find([c.pvsystem.kVA] < 500);
pv = c.pvsystem(pv_id);
%% duplicate all existing system by 5 more times
pv5 = repmat(pv,1,5);
% rename system
for i = 1:length(pv5)
    pv5(i).Name = ['newPV' sprintf('_%03.0f',i)];
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
% bus_id = bus_id(1:length(pv5));
% save('data/f909busPV.mat','bus','bus_id');

%% replace PV buses with ramdom buses
load([normalizePath('$KLEISSLLAB24-1') 'database/gridIntegration/909/f909busPV.mat'])
pbcounter=0;
rdmpv=0;

for i = 1:length(pv5)
    if (~ismember(bus{bus_id(i)},cleanBus(transbuslist))) 
        pv5(i).bus1 = bus{bus_id(i)};
    else
        pv5(i).bus1 = '0909100C';
        pbcounter=pbcounter+1;
    end
end
disp([num2str(pbcounter) ' pb with transformers'])

%% add pv to original system
c.pvsystem = [pv0 pv5];
pv = [pv0 pv5];

save([pwd '/data/f909scenario2_pv.mat'],'pv');
%% reg settings
c.transformer(end).Name = lower(c.transformer(1).Name);
c.transformer(end).sub = 'y';
c.regcontrol = dssregcontrol;
c.regcontrol(1).transformer = lower(c.transformer(1).Name);
c.regcontrol(1).Name = c.regcontrol(1).transformer ;
c.regcontrol.winding = 2;
c.regcontrol.ptratio= 12000/120/sqrt(3);
c.regcontrol(1).vreg = 121;
c.regcontrol(1).band = 1;
c.capacitor.Name = ['c' c.capacitor.Name];
c.capcontrol.Capacitor = c.capacitor.Name;
c.capcontrol.ONsetting=118;
c.capcontrol.OFFsetting=120;
c=rmfield(c,'capcontrol');
for l=1:length(c.capacitor)
    c.capacitor(l).Numsteps= 1;
end
% Energymeter.
c.energymeter = dssenergymeter;
c.energymeter.Name='substation';
c.energymeter.element=['transformer.' c.transformer(end).Name];
save('data/f909_test.mat','c')

%% add PV profils
for iii=1:3
% iii=3
    circuitPath = [pwd '/data/f909_test.mat']; % file generated upper in the page
    pvForecastPath = [normalizePath('$KLEISSLLAB24-1') 'database/gridIntegration/909/' Days(iii).day]; 
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
        summprod = summprod + c.generator(i).Kw*c.loadshape(find(ismember({c.loadshape.Name},c.generator(i).Daily))).Mult(1529);
%         prof1500(i-1) = c.loadshape(find(ismember({c.loadshape.Name},c.generator(i).Daily))).Mult(1500);
%         yop(i-1) = isequal(c.loadshape(find(ismember({c.loadshape.Name},c.generator(i).Daily))).Mult, zeros(1,2880)+0.0001);
    end
    summprod
    
    % change capacitance
%     for u=1:length(c.linecode)
%         c.linecode(u).C0=1.6*0.8;
%         c.linecode(u).C1=3.4*0.8;
%         if c.linecode(u).X1>c.linecode(u).R1
%             c.linecode(u).X0=c.linecode(u).R0;
%             c.linecode(u).X1=c.linecode(u).R1;
%         end
%     end
    % unique([c.line.C0])
    % unique([c.line.C1])
    %% Simulation no pv
    for k=1:length(c.generator)
        c.generator(k).Enabled='false';
    end
    c_bk=c;
    c.dataNoPV = dssSimulation(c,'daily',30);
    %% 15% pen
    c.generator = c_bk.generator;
    factor = 0.261;%0.5172
    for k=1:19
        c.generator(k).Enabled='true';
        c.generator(k).Kw=c.generator(k).Kw*factor;
    end
    disp(['PV penetration ' num2str(sum([c.generator{1:19}.Kw])/sum([c.load.Kw]))])
    c.dataPV15 = dssSimulation(c,'daily',30);

    sunset = 17; % 5pm 
%     c2pnames= fieldnames(c.dataPV15.Cap2Plot);
%     for k=1:length(c.dataPV15.Cap2Plot)
%         c.dataPV15.Cap2Plot.(c2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Cap2Plot.(c2pnames{k})(2880*sunset/24:end);
%     end
    t2pnames= fieldnames(c.dataPV15.Tap2Plot);
    for k=1:length(c.dataPV15.Tap2Plot)
        c.dataPV15.Tap2Plot.(t2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Tap2Plot.(t2pnames{k})(2880*sunset/24:end);
    end
    c.dataPV15.LineLoss(2880*sunset/24:end,:)=c.dataNoPV.LineLoss(2880*sunset/24:end,:);
    c.dataPV15.TotalLoss(2880*sunset/24:end,:)=c.dataNoPV.TotalLoss(2880*sunset/24:end,:);
    c.dataPV15.TotalPower(2880*sunset/24:end,:)=c.dataNoPV.TotalPower(2880*sunset/24:end,:);
    c.dataPV15.VoltMaxMin(2880*sunset/24:end,:)=c.dataNoPV.VoltMaxMin(2880*sunset/24:end,:);
    c.dataPV15.Voltage(2880*sunset/24:end,ismember(c.dataPV15.nodeName,c.dataNoPV.nodeName)) = ...
        c.dataPV15.Voltage(2880*sunset/24:end,ismember(c.dataPV15.nodeName,c.dataNoPV.nodeName));
    c.dataPV15.Capcontrol(2880*sunset/24:end,:) = c.dataNoPV.Capcontrol(2880*sunset/24:end,:);
    c.dataPV15.Regulation(2880*sunset/24:end,:) = c.dataNoPV.Regulation(2880*sunset/24:end,:);
    
    %% 30% pen
    c.generator = c_bk.generator;
factor = 0.521;
    for k=1:19
        c.generator(k).Enabled='true';
        c.generator(k).Kw=c.generator(k).Kw*factor;
    end
    disp(['PV penetration ' num2str(sum([c.generator{1:19}.Kw])/sum([c.load.Kw]))])
    c.dataPV30 = dssSimulation(c,'daily',30);

    sunset = 17; % 5pm 
%     c2pnames= fieldnames(c.dataPV30.Cap2Plot);
%     for k=1:length(c.dataPV30.Cap2Plot)
%         c.dataPV30.Cap2Plot.(c2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Cap2Plot.(c2pnames{k})(2880*sunset/24:end);
%     end
    t2pnames= fieldnames(c.dataPV30.Tap2Plot);
    for k=1:length(c.dataPV30.Tap2Plot)
        c.dataPV30.Tap2Plot.(t2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Tap2Plot.(t2pnames{k})(2880*sunset/24:end);
    end
    c.dataPV30.LineLoss(2880*sunset/24:end,:)=c.dataNoPV.LineLoss(2880*sunset/24:end,:);
    c.dataPV30.TotalLoss(2880*sunset/24:end,:)=c.dataNoPV.TotalLoss(2880*sunset/24:end,:);
    c.dataPV30.TotalPower(2880*sunset/24:end,:)=c.dataNoPV.TotalPower(2880*sunset/24:end,:);
    c.dataPV30.VoltMaxMin(2880*sunset/24:end,:)=c.dataNoPV.VoltMaxMin(2880*sunset/24:end,:);
    c.dataPV30.Voltage(2880*sunset/24:end,ismember(c.dataPV30.nodeName,c.dataNoPV.nodeName)) = ...
        c.dataPV30.Voltage(2880*sunset/24:end,ismember(c.dataPV30.nodeName,c.dataNoPV.nodeName));
    c.dataPV30.Capcontrol(2880*sunset/24:end,:) = c.dataNoPV.Capcontrol(2880*sunset/24:end,:);
    c.dataPV30.Regulation(2880*sunset/24:end,:) = c.dataNoPV.Regulation(2880*sunset/24:end,:);
    
    
    %% Real PV pen 57.6%
    c.generator = c_bk.generator;
    disp(['PV penetration ' num2str(sum([c.generator{1:19}.Kw])/sum([c.load.Kw]))])
    for k=1:19
        c.generator(k).Enabled='true';
        c.generator(k).Kw=c.generator(k).Kw;
    end
    c.dataPV58 = dssSimulation(c,'daily',30);

    sunset = 17; % 5pm 
%     c2pnames= fieldnames(c.dataPV58.Cap2Plot);
%     for k=1:length(c.dataPV58.Cap2Plot)
%         c.dataPV58.Cap2Plot.(c2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Cap2Plot.(c2pnames{k})(2880*sunset/24:end);
%     end
    t2pnames= fieldnames(c.dataPV58.Tap2Plot);
    for k=1:length(c.dataPV58.Tap2Plot)
        c.dataPV58.Tap2Plot.(t2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Tap2Plot.(t2pnames{k})(2880*sunset/24:end);
    end
    c.dataPV58.LineLoss(2880*sunset/24:end,:)=c.dataNoPV.LineLoss(2880*sunset/24:end,:);
    c.dataPV58.TotalLoss(2880*sunset/24:end,:)=c.dataNoPV.TotalLoss(2880*sunset/24:end,:);
    c.dataPV58.TotalPower(2880*sunset/24:end,:)=c.dataNoPV.TotalPower(2880*sunset/24:end,:);
    c.dataPV58.VoltMaxMin(2880*sunset/24:end,:)=c.dataNoPV.VoltMaxMin(2880*sunset/24:end,:);
    c.dataPV58.Voltage(2880*sunset/24:end,ismember(c.dataPV58.nodeName,c.dataNoPV.nodeName)) = ...
        c.dataPV58.Voltage(2880*sunset/24:end,ismember(c.dataPV58.nodeName,c.dataNoPV.nodeName));
    c.dataPV58.Capcontrol(2880*sunset/24:end,:) = c.dataNoPV.Capcontrol(2880*sunset/24:end,:);
    c.dataPV58.Regulation(2880*sunset/24:end,:) = c.dataNoPV.Regulation(2880*sunset/24:end,:);

    %% 74.4% PV pen 
    disp(['PV penetration ' num2str(sum([c.generator.Kw])/sum([c.load.Kw]))])
    for k=1:length(c.generator)
        c.generator(k).Enabled='true';
    %     c.generator(k).Kw=c.generator(k).Kw;
    end
    c.dataPV75 = dssSimulation(c,'daily',30);

    sunset = 17; % 5pm 
%     c2pnames= fieldnames(c.dataPV75.Cap2Plot);
%     for k=1:length(c.dataPV75.Cap2Plot)
%         c.dataPV75.Cap2Plot.(c2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Cap2Plot.(c2pnames{k})(2880*sunset/24:end);
%     end
    t2pnames= fieldnames(c.dataPV75.Tap2Plot);
    for k=1:length(c.dataPV75.Tap2Plot)
        c.dataPV75.Tap2Plot.(t2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Tap2Plot.(t2pnames{k})(2880*sunset/24:end);
    end
    c.dataPV75.LineLoss(2880*sunset/24:end,:)=c.dataNoPV.LineLoss(2880*sunset/24:end,:);
    c.dataPV75.TotalLoss(2880*sunset/24:end,:)=c.dataNoPV.TotalLoss(2880*sunset/24:end,:);
    c.dataPV75.TotalPower(2880*sunset/24:end,:)=c.dataNoPV.TotalPower(2880*sunset/24:end,:);
    c.dataPV75.VoltMaxMin(2880*sunset/24:end,:)=c.dataNoPV.VoltMaxMin(2880*sunset/24:end,:);
    c.dataPV75.Voltage(2880*sunset/24:end,ismember(c.dataPV75.nodeName,c.dataNoPV.nodeName)) = ...
        c.dataPV75.Voltage(2880*sunset/24:end,ismember(c.dataPV75.nodeName,c.dataNoPV.nodeName));
    c.dataPV75.Capcontrol(2880*sunset/24:end,:) = c.dataNoPV.Capcontrol(2880*sunset/24:end,:);
    c.dataPV75.Regulation(2880*sunset/24:end,:) = c.dataNoPV.Regulation(2880*sunset/24:end,:);

    %% 100% PV pen 
    c.generator = c_bk.generator;
    for k=[2:10,12:length(c.generator)]
        c.generator(k).Enabled='true';
        c.generator(k).Kw=c.generator(k).Kw*2.275;
    end
    c.generator(1).Enabled='true';
    c.generator(11).Enabled='true';
    disp(['PV penetration ' num2str(sum([c.generator{ismember({c.generator.Enabled},'true')}.Kw])/sum([c.load.Kw]))])
    c.dataPV100 = dssSimulation(c,'daily',30);

    sunset = 17; % 5pm 
%     c2pnames= fieldnames(c.dataPV100.Cap2Plot);
%     for k=1:length(c.dataPV100.Cap2Plot)
%         c.dataPV100.Cap2Plot.(c2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Cap2Plot.(c2pnames{k})(2880*sunset/24:end);
%     end
    t2pnames= fieldnames(c.dataPV100.Tap2Plot);
    for k=1:length(c.dataPV100.Tap2Plot)
        c.dataPV100.Tap2Plot.(t2pnames{k})(2880*sunset/24:end) = c.dataNoPV.Tap2Plot.(t2pnames{k})(2880*sunset/24:end);
    end
    c.dataPV100.LineLoss(2880*sunset/24:end,:)=c.dataNoPV.LineLoss(2880*sunset/24:end,:);
    c.dataPV100.TotalLoss(2880*sunset/24:end,:)=c.dataNoPV.TotalLoss(2880*sunset/24:end,:);
    c.dataPV100.TotalPower(2880*sunset/24:end,:)=c.dataNoPV.TotalPower(2880*sunset/24:end,:);
    c.dataPV100.VoltMaxMin(2880*sunset/24:end,:)=c.dataNoPV.VoltMaxMin(2880*sunset/24:end,:);
    c.dataPV100.Voltage(2880*sunset/24:end,ismember(c.dataPV100.nodeName,c.dataNoPV.nodeName)) = ...
        c.dataPV100.Voltage(2880*sunset/24:end,ismember(c.dataPV100.nodeName,c.dataNoPV.nodeName));
    c.dataPV100.Capcontrol(2880*sunset/24:end,:) = c.dataNoPV.Capcontrol(2880*sunset/24:end,:);
    c.dataPV100.Regulation(2880*sunset/24:end,:) = c.dataNoPV.Regulation(2880*sunset/24:end,:);
%% save and plot results 
    save(['data/f909_SimResults_' Days(iii).weather '.mat'],'c');
    for h=[15,30, 58,75,100]
        mkdir([pwd '/tmp/f909' num2str(h) '/'])
        result_dir=[pwd '/tmp/f909' num2str(h) '/'];
        DailyPlot(c.(['dataPV' num2str(h)]), c.dataNoPV, result_dir,...
            ['F909_' num2str(h) 'pen_' Days(iii).weather],1,[13 13.5 14]);
    end
end
%% Plot comparison between all cases
clear Data;
Data ={};
PVpen= [0, 15, 30, 58, 75, 100];
for w=1:3   %weather 
    resu_file=['data/f909_SimResults_' Days(w).weather '.mat'];
    c=load(resu_file);c=c.c;
    Data{w,1} = c.dataNoPV;
    for g=2:length(PVpen) 
        Data{w,g} = c.(['dataPV' num2str(PVpen(g))]);
    end
end


leg={'cloudy','overcast','clear'};
outPath = 'tmp/f909';
dataplotPVpen( PVpen, Data ,leg, outPath)