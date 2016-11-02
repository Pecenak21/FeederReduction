%% summary
% 1. add all 432 PV systems (including 2 big) - save as data/f520_432pv.mat
% 2. add pv profiles (disaggregated + aggregated using 1st profile)
% 3. fixed 30 empty pv profiles (resulted from overlapping footprints) -> no longer needed if forecast is rerun 
%       This file ('$KLEISSLLAB24-1') 'database\gridIntegration\USI_1_2/PVchanges.mat) contains the overlapping sites
%		Comment out this when forecast is rerun (4 line from the filename above)
% 4. created generator + rm pvsystem
% 5. VR settings
% 6. generate 'data/f520Solved.mat' files including the results of 5 cases with diff penetrations (always with big PV systems) 
% 7. run simulation for each case with 3 cloud conditions
% 8. plot comparison figures
% 
% config: run days '20121214', '20121218','20121219' for forecast. 
%			use 1 day '6/2/2013' to get loadshape from SDGE Point Loma
%		ISSUE: simulation only run for '6/2/2013' to get loadshape from SDGE Point Loma and no other days 

%% Presetting
clear;
c=load('data/f520.mat','c');
c=c.c;
Days = struct('daysForForecast',{'20121214', '20121218','20121219' },'weather',{'cloudy','overcast','clear'},...
    'dateForLoadshapeFromSDGE',{'12/19/2012','12/19/2012','12/19/2012'}); %{'12/14/2012','12/18/2012','12/19/2012'}; %{'5/25/2013','4/26/2013','6/2/2013'};% %
c.pvsystem = [];
pv=load([normalizePath('$KLEISSLLAB24-1') 'database\gridIntegration\Fallbrook setting\Fallbrook_Scenario2\scenario2_pvProfiles.mat']); pv=pv.pv;
c.pvsystem=dsspvsystem;
for i=1:length(pv)
    c.pvsystem(i) = dsspvsystem('Name',['pv' num2str(i)]);
    c.pvsystem(i).bus1 = pv(i).bus1;
    c.pvsystem(i).kVA = pv(i).kVA;
    c.pvsystem(i).Pmpp = pv(i).kVA;
end
%     %% Change regcontrolers' name
for tra = 1:length(c.transformer)
    asd = regexp(c.transformer(tra).Name,'\_','split');
    c.transformer(tra).Name = ['t' lower([asd{:}])];
    c.regcontrol(tra).Name = c.transformer(tra).Name;
    c.regcontrol(tra).transformer = c.transformer(tra).Name;
end
% Energymeter.
c.energymeter = dssenergymeter;
c.energymeter.Name='substation';
c.energymeter.element=['transformer.' c.transformer(7).Name];
save('data/f520_432pv.mat','c')    
% %% Loop
for iii=[1,2,3]
%     %% Add pv profil
%     iii=3;
    circuitPath = [pwd '/data/f520_432pv.mat']; % file generated upper in the page
    pvForecastPath = [normalizePath('$KLEISSLLAB24-1') 'database/gridIntegration/USI_1_2/' Days(iii).daysForForecast]; 
    desag_ag = [1 1]; 
    if (desag_ag==0); profil_type = 'aggregated'; else profil_type = 'disaggregated'; end
    outPath = ['tmp/'];
    addBigPlant=[0 0];
    plotOption = 0;

    fprintf(['Creating ' profil_type ' pv profil:\n'])
    [c p] = AddPVForecast(circuitPath, pvForecastPath, desag_ag, outPath, addBigPlant, plotOption);
%     %% Fix forecast problems
    pvchanges = load([normalizePath('$KLEISSLLAB24-1') 'database\gridIntegration\USI_1_2/PVchanges.mat']);pvchanges=pvchanges.PVchanges;
    for w=1:length(pvchanges)
        c.pvsystem(pvchanges(w,1)).daily = c.pvsystem(pvchanges(w,2)).daily;
    end
%     for w=1:length(pvchanges)
%         c.loadshape(ismember({c.loadshape.Name},c.pvsystem(pvchanges(length(pvchanges)+1-w,1)).daily)) =[];
%     end
%     %% PVGene
    offset=length(c.generator);
    for i = 1:length(c.pvsystem)
        % Panel kW = Pmpp (in kW @1kW/m2 and 25 C) * Irradiance (in kW/m2) * Factor(@actual T)
        c.generator(i+offset).Name = ['gen_' num2str(i)];
        c.generator(i+offset).Kw = c.pvsystem(i).kVA;
        c.generator(i+offset).Kvar = 0;
        c.generator(i+offset).kv = 12;
        c.generator(i+offset).phases = 1;
        c.generator(i+offset).Vmaxpu = 1.6;
        c.generator(i+offset).Vminpu = 0.8;
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
    c.generator(1).Enabled = 'false';
    % for k=46:length(c.generator)
    %     c.generator(k).Enabled='false';
    % end
%     %% Settings
    for reg=1:length(c.regcontrol)
        c.regcontrol(reg).winding = 2;
    end
    c.regcontrol(1).vreg = 120.5;
    c.regcontrol(1).band = 1.5;
    c.regcontrol(2).vreg = 120.5;
    c.regcontrol(2).band = 1.5;
    c.regcontrol(3).vreg = 120.5;
    c.regcontrol(3).band = 1.5;
    c.regcontrol(4).vreg = 120.5;
    c.regcontrol(4).band = 1.5;
    c.regcontrol(5).vreg = 120.5;
    c.regcontrol(5).band = 1.5;
    c.regcontrol(6).vreg = 120.5;
    c.regcontrol(6).band = 1.5;
    c.regcontrol(7).vreg = 121;
    c.regcontrol(7).band = 1;
%     c.regcontrol(2)=[];
    % c.transformer(4).MaxTap = 1.1;
    % c.transformer(5).MaxTap = 1.1;
    c.capcontrol(3).Name = ['cc_' cleanBus(c.capacitor(3).Name)];
    c.capcontrol(3).Capacitor = c.capacitor(3).Name;
    c.capcontrol(3).Element = 'line.05201643_05201643A';
    c.capcontrol(3).Type ='voltage';
    c.capcontrol(3).PTRatio = 57.75;
    c.capcontrol(3).Vmax = 126; c.capcontrol(3).Vmin = 117;
    

    for cap=1:length(c.capcontrol)
        c.capcontrol(cap).OFFsetting = 121;
        c.capcontrol(cap).ONsetting = 119;
        c.capcontrol(cap).CTRatio = 300;
    end
    c.capcontrol(3).OFFsetting = 120;
    c.capcontrol(3).ONsetting = 118;

%     %% Add loadshape
    load([pwd '/data\load.mat'])
%     %% Create load profiles 
    [ c,  p] = AddLS( c , pwd);%, date{iii} );
% % %     %aggregated loadshapes
% % %     l1 = la;
% % %     [x idx] = ismember({l1.Name}',strcat('Load.',{c.load.Name}'));
% % %     tic
% % %     dt = 30/3600/24;
% % %     t = dt:dt:1;
% % %     if isfield(c,'loadshape');offset = length(c.loadshape);else offset=0;end
% % %     for i = 1:length(l1)
% % %         % get id of loadshape considering the offset
% % %         id = offset + i;
% % %         c.loadshape(id) = dssloadshape;
% % %         c.loadshape(id).Name = ['loadshape_' c.load(idx(i)).Name];
% % %         c.loadshape(id).sInterval = 30;
% % %         c.loadshape(id).Npts = length(t);
% % %         % set up multiplier; assume the load profile starts at 12am in the morning
% % %         m = [l1(i).x12_00_AM l1(i).x1_00_AM l1(i).x2_00_AM l1(i).x3_00_AM l1(i).x4_00_AM l1(i).x5_00_AM...
% % %             l1(i).x6_00_AM l1(i).x7_00_AM l1(i).x8_00_AM l1(i).x9_00_AM l1(i).x10_00_AM l1(i).x11_00_AM l1(i).x12_00_N...
% % %             l1(i).x1_00_PM l1(i).x2_00_PM l1(i).x3_00_PM l1(i).x4_00_PM l1(i).x5_00_PM l1(i).x6_00_PM...
% % %             l1(i).x7_00_PM l1(i).x8_00_PM l1(i).x9_00_PM l1(i).x10_00_PM l1(i).x11_00_PM l1(i).x12_00_AM];
% % %         % interpolate to get 30 second interval data
% % %         c.loadshape(id).Mult = interp1(0:1/24:1,m,t,'cubic');
% % %         c.load(idx(i)).Daily = c.loadshape(id).Name;
% % %         if ~(i==1)
% % %             totalLS = totalLS + c.loadshape(id).Mult*c.load(i).kVA;
% % %         else 
% % %             totalLS = c.loadshape(id).Mult*c.load(i).kVA;
% % %         end
% % %     end
    % figure, plot(0:24/2880:24-24/2880,totalLS);
% c.loadshape(end).Mult = c.loadshape(end).Mult*0.8;%zeros(1,2880);%c.loadshape(1).Mult*0.9;
	c=rmfield(c,'capcontrol');
	for l=1:length(c.capacitor)
		c.capacitor(l).Numsteps= 1;
	end
    save([pwd '/data/f520_test.mat'],'c')


%     %% Plot
%     sumpv=0;
%     for k=2:length(c.generator)
% %         if isequal(c.generator(k).Enabled,'true')
%             sumpv= sumpv + c.generator(k).Kw*c.loadshape(find(ismember({c.loadshape.Name},c.generator(k).Daily))).Mult;
% %         end
%     end
%     figure, plot(sumpv/1000)

%     %% Launch Simulation
    % c=addMonitors(c);
    % circuitVisualizerTEST(c)

%     %% Cases (load max rating : 11.123MW, 13.192MVA) 
    %kvaload = sum(sqrt([c.load.Kw].*[c.load.Kw]+[c.load.Kvar].*[c.load.Kvar]))
%     %% Case 1: Normal 45 Pv systems - total pv rating: 2.2954MW - PV PENETRATION: 17.40%
    c= load('data/f520_test');c=c.c;
    cbk=c;
    for k=1:length(c.generator)
        c.generator(k).Enabled='false';
    end
    disp('No PV:')
    c.dataNoPV =dssSimulation(c,'daily',30);
    
    for k=2:46
        c.generator(k).Enabled='true';
    end
    for k=47:length(c.generator)
        c.generator(k).Enabled='false';
    end
    disp('17% PV pen.:')
    c.dataPV17 =dssSimulation(c,'daily',30);
    %
    for k=2:201
        c.generator(k).Enabled='true';
    end
    for k=202:length(c.generator)
        c.generator(k).Enabled='false';
    end
    disp('30% PV pen.:')
    c.dataPV30 =dssSimulation(c,'daily',30);
    
    sFactor=1.3;
    for k=2:length(c.generator)
        c.generator(k).Enabled='true';
    end
    for k=2:44
        c.generator(k).Kw=c.generator(k).Kw*sFactor;
    end
    for k=47:length(c.generator)
        c.generator(k).Kw=c.generator(k).Kw*sFactor;
    end
    disp('51% PV pen.:')
    c.dataPV51 =dssSimulation(c,'daily',30);
    
    c.generator=cbk.generator;
    sFactor=2.15;
    for k=2:length(c.generator)
        c.generator(k).Enabled='true';
    end
    for k=2:44
        c.generator(k).Kw=c.generator(k).Kw*sFactor;
    end
    for k=47:length(c.generator)
        c.generator(k).Kw=c.generator(k).Kw*sFactor;
    end
    disp('75% PV pen.:')
    c.dataPV75 =dssSimulation(c,'daily',30);
    
    c.generator=cbk.generator;
    sFactor=3.1;
    for k=2:length(c.generator)
        c.generator(k).Enabled='true';
    end
    for k=2:44
        c.generator(k).Kw=c.generator(k).Kw*sFactor;
    end
    for k=47:length(c.generator)
        c.generator(k).Kw=c.generator(k).Kw*sFactor;
    end
    disp('100% PV pen.:')
    c.dataPV100 =dssSimulation(c,'daily',30);
%     save('data/f520Solved','c');
    
    namelist={'dataPV17','dataPV30','dataPV51','dataPV75','dataPV100'};
    
    for l=1:5
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
        mkdir([pwd '/tmp/f520SimCase' namelist{l}(7:end) '/'])
        result_dir=[pwd '/tmp/f520SimCase' namelist{l}(7:end) '/'];
        DailyPlot(c.(namelist{l}), c.dataNoPV, result_dir,...
            ['F520_Case' namelist{l}(7:8) '_' Days(iii).weather],1,[13 13.5 14])
    end
%     %% Case 2: 200 Pv systems - total pv rating: 3.035MW - PV PENETRATION: 30.6%
   
%     %% Case 3: 432 PV systems - total pv rating: 6.756MW - PV PENETRATION: 51.2%
   
%      %% Case 4: 432 PV systems - total pv rating: 8.3946MW - PV PENETRATION: 75.5%
    
%      %% Case 5: 432 PV systems - total pv rating: 11.221MW - PV PENETRATION: 100.9%

%     %% Testing
% %     
% %     sumpv=0;
% %     scale_factor=3.1;
% %     for k=2:44%length(c.generator)
% % %         if isequal(c.generator(k).Enabled,'true')
% %             sumpv=sumpv+c.generator(k).Kw*scale_factor;
% % %         end
% %     end
% %     for k=47:length(c.generator)
% % %         if isequal(c.generator(k).Enabled,'true')
% %             sumpv=sumpv+c.generator(k).Kw*scale_factor;
% % %         end
% %     end
% %     sumpv = sumpv+c.generator(45).Kw+c.generator(46).Kw
% %     sumpv/11122
% %     
% %     %% factor 1.3 -> 52.7 
% %     %% 200PV -> 30.18 
% %     
% %     %% factor 2.15 -> 8.3946 75.5
% %     %% factor 3.1 -> 11.221 100.9
% %     
	save([pwd '/data/f520_SimResults_' Days(iii).weather '.mat'],'c');
end
% % % for p=2:length(c.generator)
% % %     blbl(p)=ismember(c.generator(p).Daily,{c.loadshape.Name});
% % % end

%% Plot comparison between all cases
clear Data;

PVpen= [0, 17, 30, 51, 75, 100];
Data ={};
for w=1:3   %weather 
    c=load([pwd '/data/f520_SimResults_' Days(w).weather '.mat']);c=c.c;
    for g=1:5 % number of cases
        if g==1,
            Data{w,g} = c.dataNoPV;
        end
        Data{w,g+1} = c.(['dataPV' num2str(PVpen(g+1))]);
    end
end


leg={'cloudy','overcast','clear'};
outPath = 'tmp/f520';
dataplotPVpen( PVpen, Data ,leg, outPath)

