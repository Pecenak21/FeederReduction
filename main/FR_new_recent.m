%New FR test code
clear all
Full=tic;
close all

useSaved=1;
%% Todo
%% Reading in
%Make it so the code reads the .dss file, if it has read the same file
%since no changes, load in variables (Ybus, PV, Load, etcc...)
%otherwise, grind through it all.
%Do this by date modified in DIr, we will probably have to save the last
%tiem we used this .dss file somewhere

%% Load vminpu
%Right nwo there is no way to account for different load Vminpu in
%aggregation.

%% capacitor subtraction from Ybus
%Just like transformer, we need to subtract our capascitors that we are
%keepign in the system from the ybus so that they are not written as
%reactors. we need to remove from Ybus in case any shutn capacitance is
%written into Ybus.

%% VR
%Right now we keep all VR, maybe we can remove any that don't have buses
%downstream. We probably can actually, use same code as one used on trf

p='C:\Program Files\OpenDSS\EPRITestCircuits\ckt24\master_ckt24.dss';
	if ~exist(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '.mat']) || useSaved
		circuit=dssparse(p);
		save(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '.mat'],'circuit')
	else
		load(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '.mat'])
	end
%% start code
if regexp(feeder,'ieee')
	tic
	p='C:\Program Files\OpenDSS\EPRITestCircuits\ckt24\master_ckt24.dss';
% 	p='c:\users\zactus\FeederReductionRepo\8500-Node\Master.dss';
% 			p='C:\Users\Zactus\Documents\OpenDSS\OpenDSS\IEEETestCases\13Bus\IEEE13Nodeckt_zack.dss';
% 				p='c:\users\zactus\FeederReductionRepo\123Bus\IEEE123Master.dss';
	% 		p='C:\Users\Zactus\Documents\OpenDSS\OpenDSS\IEEETestCases\34Bus\ieee34_1.dss';
	% 	circuit=dssparse('C:\Users\Zactus\Documents\OpenDSS\OpenDSS\IEEETestCases\13Bus\IEEE13Nodeckt_zack.dss');
	if ~exist(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '.mat']) || useSaved
		circuit=dssparse(p);
		save(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '.mat'],'circuit')
	else
		load(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '.mat'])
	end
	toc
	
	orig=tic;
	if ~exist(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_sim.mat'])  || useSaved
		[ data1 ] = dssSimulation( circuit,[],1,[],[],0);
		FullTime=toc(orig);
		save(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_sim.mat'],'data1')
	else
		load(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_sim.mat'])
	end
	
	o = actxserver('OpendssEngine.dss');
	dssText = o.Text; dssText.Command = 'Clear'; cDir = pwd;
	dssText.Command = ['Compile "' p '"'];
	dssCircuit = o.ActiveCircuit;
	circuit.buslist.id=regexprep(dssCircuit.AllBUSNames,'-','_');
	circuit.buslist.coord=zeros(length(circuit.buslist.id),2);
	delete(o);
	
	
elseif strcmp(feeder,'Alpine')
	circuit=load('c:\users\zactus\gridIntegration\results\Alpine_s1b_refined.mat');
	circuit=circuit.c;
	
	%error
	for j=1:length(circuit.pvsystem)
		if length(regexp(circuit.pvsystem(j).bus1,'\.','match'))==1 %single phase
			circuit.pvsystem(j).kV=12/sqrt(3);
		else
			circuit.pvsystem(j).kV=12;
		end
	end
	% 	circuit=rmfield(circuit,'load');
	% 	circuit=rmfield(circuit,'pvsystem');
	
	[ data1 ] = dssSimulation( circuit,[],1,[],[],0);
	% 	criticalBuses=circuit.buslist.id(1:end-3);
	criticalBuses=circuit.buslist.id([11,13,14,456,603]);
	
	
end
% criticalBuses=circuit.buslist.id([1,100,1000,2000,3000:end]);
% criticalBuses=circuit.buslist.id([1,100]);

% criticalBuses=circuit.buslist.id;
% criticalBuses=buslist;

% criticalBuses=regexprep(criticalBuses,'-','_');

% criticalBuses=circuit.buslist.id(round((length(circuit.buslist.id)-2)*rand(kk,1)+1));
% criticalBuses=unique(criticalBuses);
% while length(criticalBuses)<kk
% 	cnode_tmp=circuit.buslist.id(round((length(circuit.buslist.id)-2)*rand(kk,1)+1));
% 	criticalBuses=[criticalBuses;cnode_tmp];
% 	criticalBuses=unique(criticalBuses);
% end

% criticalBuses=circuit.buslist.id(1:end-3);
% criticalBuses=circuit.buslist.id;


% criticalBuses=circuit.buslist.id([11,13,14,456,603]);
% criticalBuses=circuit.buslist.id(find(ismemberi(circuit.buslist.id,{'633','634','671'})));
% criticalBuses=circuit.buslist.id(find(ismemberi(circuit.buslist.id,{'632'})));


%% check which variables we have and set flags
fprintf('\n-------------------------------------------------------------\n')
fprintf('\n                     Initializing circuit        \n')
fprintf('\n-------------------------------------------------------------\n')

[Flag]=isfield(circuit,{'load','pvsystem','capacitor','transformer','regcontrol','capcontrol','reactor'});
circuit_orig=circuit;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Y-BUS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
fprintf('\nGetting Ybus: ')
if ~exist(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Ybus.mat'])  || useSaved
	%remove shunt terms from load/pv
	if Flag(2)
		circuitBase=rmfield(circuit,'pvsystem');
	else
		circuitBase=circuit;
	end
	if Flag(1)
		circuitBase=rmfield(circuitBase,'load');
	end
	if Flag(5)
		circuitBase=rmfield(circuitBase,'regcontrol');
	end
	if Flag(4)
		circuit_noXfrmr=rmfield(circuitBase,'transformer');
	end
	
	[~, ~, Ycomb_noXfrmr, Ybus_noXfrmr, ~,~]=getYbus(circuit_noXfrmr);
	[YbusOrderVect, YbusPhaseVect, Ycomb, Ybus, buslist,dssCircuit]=getYbus(circuitBase);
	
	
	missingNode=Ycomb(find(~ismember(Ycomb,Ycomb_noXfrmr)));
	%Add in any missign buses to new Ybus
	Ycomb_noXfrmr=[Ycomb_noXfrmr;missingNode];
	Ybus_noXfrmr=[Ybus_noXfrmr zeros(length(Ybus_noXfrmr),length(missingNode)); zeros(length(missingNode),length(Ybus_noXfrmr)+length(missingNode))];
	
	%Make sure Order of Ybus matches order of Ybus
	%This function takes keySet and assigns an order to that (i.e. 'A'=1,
	%'B'=2...
	%Then it tells you the order of valuSet in terms of keySet. (i.e Valueset
	%=['B', 'A'], values(mapObj,valueSet)= [2 1];
	valueSet = lower(Ycomb_noXfrmr);
	keySet = lower(Ycomb);
	mapObj = containers.Map(keySet,1:length(keySet));
	[~,Yorder]=sort(cell2mat(values(mapObj,valueSet)));
	Ycomb_noXfrmr=Ycomb_noXfrmr(Yorder);
	Ybus_noXfrmr=Ybus_noXfrmr(Yorder,Yorder);
	
	YofXfrmrs=Ybus-Ybus_noXfrmr;
	
	circuit_orig.Ybus=Ybus;
	Zbus=inv(Ybus);
	circuit_orig.Zbus=Zbus;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% %Get PU ybus
	baseKv=unique(cell2mat([circuit.transformer{:}.kV]));
	RealV=dssCircuit.AllBusVmag*sqrt(3)/1000;
	if length(baseKv)>1
		baseKvMat=repmat(baseKv,length(RealV),1);
		volt_base=[];
		VoltDiff = bsxfun(@minus,baseKvMat,RealV');
		[~,Ind]=min(abs(VoltDiff),[],2);
		volt_base=baseKvMat(sub2ind(size(baseKvMat),[1:size(Ind,1)]',Ind));
	else
		volt_base=repmat(baseKv,length(RealV),1);
	end
	
	%Make sure Order of voltages matches order of Ybus
	valueSet = lower(Ycomb);
	keySet = lower(dssCircuit.AllnodeNames);
	mapObj = containers.Map(keySet,1:length(keySet));
	volt_base=volt_base(cell2mat(values(mapObj,valueSet)));
	
	power_base=1;
	Zbase=volt_base*volt_base'/power_base;
	Ybase=(1./Zbase);
	
	
	
	%get voltage vector
	% 	V3=volt_base./sqrt(3);
	% 	V3(find(ismemberi(YbusPhaseVect,2)))=V3(find(ismemberi(YbusPhaseVect,2)))*exp(-1i*2*pi/3);
	% 	V3(find(ismemberi(YbusPhaseVect,3)))=V3(find(ismemberi(YbusPhaseVect,3)))*exp(1i*2*pi/3);
	%account for phase angle in voltage
	VOLT=dssCircuit.AllBusVolts;
	ineven=2:2:length(VOLT); inodd=1:2:length(VOLT);
	VOLT=VOLT(inodd)+1i*VOLT(ineven);
	VOLT=VOLT(cell2mat(values(mapObj,valueSet)));
	ANGLE=round((angle(VOLT))/(pi/6))*(pi/6);
	
	
	V2=volt_base./sqrt(3).*exp(1i*ANGLE');
	V2=diag(V2);
	save(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Ybus.mat'],'Ybus','YbusOrderVect','YbusPhaseVect','Ycomb','buslist','YofXfrmrs','Zbus','circuit_orig','V2');
else
	load(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Ybus.mat']);
end
t_=toc;
fprintf('time elapsed %f\n',t_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Buslist    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
fprintf('\nGetting Buslist: ')
%get rid of this using key vals!!!!!!!!!!!!!!!
%for IEEE8500, since buses are missing from buslist.id
criticalBuses=buslist(find(~cellfun(@isempty,regexp(buslist,'sx'))));

criticalBuses=regexprep(criticalBuses,'-','_');
for ii=1:length(buslist)
	Ind=find(strcmpi(buslist(ii),YbusOrderVect))';
	Node_number(Ind)=ii;
end
% for j=1:length(criticalBuses)
% 	criticalNumbers(j)=find(strcmpi(buslist,lower(criticalBuses(j))));
% end
criticalNumbers=find(ismemberi(buslist,criticalBuses));
criticalNumbers=unique(criticalNumbers);
criticalBuses=buslist(criticalNumbers);

if ~iscolumn(criticalNumbers)
	criticalNumbers=criticalNumbers';
end

%Add sourcebus
t_=toc;
fprintf('time elapsed %f\n',t_)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     TOPO    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
fprintf('\nGetting Topology: ')
if ~exist(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_initTopo.mat'])  || useSaved
	topo=zeros(max(Node_number),4);
	generation{1,1}=[];clear generation;  generation{1,1}=1; generation{1,4}=[];generation{1,5}=0;
	parent=1;
	topo(parent,1)=parent;
	[topo,generation]=topology_detect_large(topo,generation,Ybus,parent,Node_number);
	topo_view=topo;
	topo_view(find(topo_view(:,1)==0)',:)=[];
	c_new=0;
	save(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_initTopo.mat'],'generation','topo');
else
	load(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_initTopo.mat']);
end
t_=toc;
fprintf('time elapsed %f\n',t_)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Capacitors    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section is to keep the capacitors in the grid by adding the nodes that they are connected to.
if Flag(3)
	tic
	fprintf('\nAdding Capacitor Nodes: ')
	if ~exist(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Cap.mat'])  || useSaved
		bus1=[];
		cap=circuit.capacitor;
		for ii=1:length(cap)
			bus1{ii}=cap(ii).bus1;
			bus2{ii}=cap(ii).bus2;
		end
		if ~iscolumn(bus1)
			bus1=bus1'; bus2=bus2';
		end
		bus1=strtok(bus1,'\.'); CapNum=find(ismemberi(buslist,bus1));
		criticalNumbers=[criticalNumbers; CapNum];
		bus2=strtok(bus2,'\.'); CapNum=find(ismemberi(buslist,bus2));
		criticalNumbers=[criticalNumbers; CapNum];
		
		criticalBuses=buslist(criticalNumbers);
		clear bus1 bus2
		
		if Flag(6)
			capcon=circuit.capcontrol;
			for ii=1:length(capcon)
				buses=regexp(capcon(ii).Element,'\.','split');
				if strcmp(buses{1},'line')
					LineNo=find(ismemberi({circuit.line{:}.Name},buses{2}));
					capConBuses1{ii}=strtok(circuit.line(LineNo).bus1,'.');
					capConBuses2{ii}=strtok(circuit.line(LineNo).bus2,'.');
				end
			end
			CapNum=reshape([find(ismemberi(buslist,capConBuses1)), find(ismemberi(buslist,capConBuses2))],[],1);
			criticalNumbers=[criticalNumbers; CapNum];
		end
		criticalBuses=buslist(criticalNumbers);
		clear bus1 bus2
		
		save(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Cap.mat'],'criticalNumbers','criticalBuses')
	else
		load(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Cap.mat'])
	end
	t_=toc;
	fprintf('time elapsed %f\n',t_)
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     LOADS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%here we preserve by load type!
if Flag(1)
	fprintf('Initialize loads: ')
	tic
	if ~exist(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Load.mat'])  || useSaved
		ld=circuit.load;
		S_LD=zeros(length(YbusOrderVect),8);
		% 	kV_LD=zeros(length(YbusOrderVect),length(YbusOrderVect));
		Ra=(1/sqrt(3))*(cosd(-30)+1i*sind(-30));
		Rb=(1/sqrt(3))*(cosd(30)+1i*sind(30));
		
		for j=1:length(ld)
			if ld(j).kVAr>0 && ld(j).kW>0
				s_ld=(ld(j).Kw+1i*ld(j).kVAr);
			elseif ld(j).kVAr>0 && ld(j).pf>0
				Q=ld(j).kVAr; pf=ld(j).pf; P=Q/sqrt(((1-pf^2))/pf^2);
				s_ld=P+1i*Q;
			elseif ld(j).kW>0 && ld(j).pf>0
				P=ld(j).kW; pf=ld(j).pf; Q=P*sqrt((1-pf^2)/pf^2);
				s_ld=P+1i*Q;
			elseif ld(j).kVA>0 && ld(j).pf>0
				S=ld(j).kVA; pf=ld(j).pf; P=pf*S; Q=sqrt(S^2-P^2);
				s_ld=P+1i*Q;
			elseif ld(j).kVA>0 && ld(j).kW>0
				S=ld(j).kVA; P=ld(j).kW; Q=sqrt(S^2-P^2);
				s_ld=P+1i*Q;
			elseif ld(j).kVA>0 && ld(j).kVAr>0
				S=ld(j).kVA; Q=ld(j).kVAr; P=sqrt(S^2-Q^2);
				s_ld=P+1i*Q;
			end
			
			
			if strcmpi(ld(j).conn,'delta')
				stop=1;
			end
			if isempty(regexp(ld(j).bus1,'\.','match')) %3 phase
				Ind=find(ismemberi(YbusOrderVect,ld(j).bus1));
				S_LD(Ind,ld(j).model)=S_LD(Ind,ld(j).model)+s_ld/3;
				% 			kV_LD(Ind,Ind)=ld(j).kV/sqrt(3);
			elseif length(regexp(ld(j).bus1,'\.','match'))>1 %2 phase
				name=regexp(ld(j).bus1,'\.','split');
				numPhases=length(name)-1;
				% 			kV_LD(Ind,Ind)=ld(j).kV/sqrt(3);
				if numPhases==2 && strcmpi(ld(j).conn,'delta')
					%Figure out ratio of Y connected loads and assign current
					%correctly
					Ind1=find(ismemberi(Ycomb,[name{1} '.' name{2}]));
					Ind2=find(ismemberi(Ycomb,[name{1} '.' name{3}]));
					
					if strcmpi([name{2} '.' name{3}],'1.2')||strcmpi([name{2} '.' name{3}],'2.3')||strcmpi([name{2} '.' name{3}],'3.1')
						S_LD(Ind1,ld(j).model)=S_LD(Ind1,ld(j).model)+Ra*s_ld;
						S_LD(Ind2,ld(j).model)=S_LD(Ind2,ld(j).model)+Rb*s_ld;
					elseif strcmpi([name{2} '.' name{3}],'2.1')||strcmpi([name{2} '.' name{3}],'3.2')||strcmpi([name{2} '.' name{3}],'1.3')
						S_LD(Ind1,ld(j).model)=S_LD(Ind1,ld(j).model)+Rb*s_ld;
						S_LD(Ind2,ld(j).model)=S_LD(Ind2,ld(j).model)+Ra*s_ld;
					end
				else
					for ii=2:length(name)
						Ind=find(ismemberi(Ycomb,[name{1} '.' name{ii}]));
						S_LD(Ind,ld(j).model)=S_LD(Ind,ld(j).model)+s_ld/numPhases;
					end
				end
			else %1 phase
				Ind=find(ismemberi(Ycomb,ld(j).bus1));
				S_LD(Ind,ld(j).model)=S_LD(Ind,ld(j).model)+s_ld;
			end
			
		end
		
		save(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Load.mat'],'S_LD')
	else
		load(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Load.mat'])
	end
	t_=toc;
	fprintf('%.2f sec\n',t_)
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PV    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Flag(2)
	fprintf('Initialize PV systems: ')
	tic
	
	if ~exist(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_PV.mat'])  || useSaved
		pv=circuit.pvsystem;
		S_PV=zeros(length(YbusOrderVect),1);
		% 	kV_PV=zeros(length(YbusOrderVect),length(YbusOrderVect));
		
		for j=1:length(pv)
			P=pv(j).pf*pv(j).kVA; Q=sqrt(pv(j).kVA^2-P^2);
			pvPower=(P+1i*Q);
			
			
			if isempty(regexp(pv(j).bus1,'\.','match')) %3 phase
				Ind=find(ismemberi(YbusOrderVect,pv(j).bus1));
				S_PV(Ind)=S_PV(Ind)+pvPower/3;
				% 			kV_PV(Ind,Ind)=pv(j).kV/sqrt(3);
			elseif length(regexp(pv(j).bus1,'\.','match'))>1 %2 phase
				% 			kV_PV(Ind,Ind)=pv(j).kV/sqrt(3);
				name=regexp(pv(j).bus1,'\.','split');
				if length(name)==3 && strcmpi(pv(j).conn,'delta')
					%Figure out ratio of Y connected loads and assign current
					%correctly
					Ind1=find(ismemberi(Ycomb,[name{1} '.' name{2}]));
					Ind2=find(ismemberi(Ycomb,[name{1} '.' name{3}]));
					
					if strcmpi([name{2} '.' name{3}],'1.2')||strcmpi([name{2} '.' name{3}],'2.3')||strcmpi([name{2} '.' name{3}],'3.1')
						S_PV(Ind1)=S_PV(Ind1)+Ra*pvPower;
						S_PV(Ind2)=S_PV(Ind2)+Rb*pvPower;
					elseif strcmpi([name{2} '.' name{3}],'2.1')||strcmpi([name{2} '.' name{3}],'3.2')||strcmpi([name{2} '.' name{3}],'1.3')
						S_PV(Ind1)=S_PV(Ind1)+Rb*pvPower;
						S_PV(Ind2)=S_PV(Ind2)+Ra*pvPower;
					end
				else
					for ii=2:length(name)
						Ind=find(ismemberi(Ycomb,[name{1} '.' name{ii}]));
						S_PV(Ind)=S_PV(Ind)+pvPower/pv(j).phases;
					end
				end
			else %1 phase
				% 			kV_PV(Ind,Ind)=pv(j).kV;
				Ind=find(ismemberi(Ycomb,pv(j).bus1));
				S_PV(Ind)=S_PV(Ind)+pvPower;
			end
		end
		% 	kV_PV(find(YbusPhaseVect==2))=kV_PV(find(YbusPhaseVect==2)).*exp(1i*120*(pi/180));
		% 	kV_PV(find(YbusPhaseVect==3))=kV_PV(find(YbusPhaseVect==3)).*exp(-1i*120*(pi/180));
		save(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_PV.mat'],'S_PV')
	else
		load(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_PV.mat'])
	end
	t_=toc;
	fprintf('%.2f sec\n',t_)
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    X-frmrs    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Handling Transformers: ')
tic
if Flag(4)

if ~exist(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Xfrmrs.mat'])  || useSaved
	
	trf=circuit.transformer;
		if Flag(5)
			rgc=circuit.regcontrol;
			for ii=1:length(rgc)
				trfName=rgc(ii).transformer;
				Trf_Ind(ii)=find(ismemberi({trf{:}.Name},trfName));
				buses=trf(Trf_Ind(ii)).buses;
				busNames=strtok(buses,'\.');
				criticalNumbers=[criticalNumbers; find(ismemberi(buslist,busNames))];
			end
		end
trfKeep=[];
	for ii=1:length(trf)
		[bus,phase]=strtok(trf(ii).Buses,'.');
		trfBus(ii,1:length(bus))=bus;
		if any(~ismember(phase(2:end),phase(1))) || any(ismember(trf(5).conn,'delta'))%Check to see if phases stay the same across Xfrmr or if delta conn
			trfKeep(ii)=ii;
		end		
	end
	% %
	% % %Get buses downstream of Vreg to allocate Vbase to load and PV
	% % trfDown=trfBus(:,2:end); trfDown=trfDown(:); trfDown=trfDown(find(~cellfun(@isempty,trfDown)));
	% % [~,trfDownInd]=ismemberi(trfDown,buslist);
	% % [trfDownInd, IA]=unique(trfDownInd); trfDown=trfDown(IA);
trfDown=trfBus(:,2:end); trfDown=trfDown(:); trfDown=trfDown(find(~cellfun(@isempty,trfDown)));
[~,trfDownInd]=ismemberi(trfDown,buslist);
trfBus(cellfun(@isempty,trfBus))={'noBus'};
	
	% trf_phase=[circuit.transformer{:}.Phases];
	% trf_phase(find(trf_phase==3))=4; trf_phase(find(trf_phase==1))=3; trf_phase(find(trf_phase==4))=1;
	% trf_kV=trf_kV(IA); %trf_kV=trf_kV.*sqrt(trf_phase);
	
	%code to make all parents into a single matrix
	c = generation(criticalNumbers,4)'; lens = sum(cellfun('length',c),1); ParentsMat = ones(max(lens),numel(lens)); ParentsMat(bsxfun(@le,[1:max(lens)]',lens)) = vertcat(c{:});
	
	%list of all xfrmrs who are candidates to be removed
	candidateToRemove=find(cellfun(@isempty,regexp(trf(:).Name,'Sub','match')));
	
	%find xfrms associatted with CB
	XfrmrKeep=unique(trfDownInd(find(ismemberi(trfDownInd,ParentsMat)))); 
	[row, ~]=ind2sub(size(trfBus),find(ismemberi(trfBus(:,2:end),buslist(XfrmrKeep)))); XfrmrKeepInd=unique(row);
% 	XfrmrKeepInd=[find(ismemberi(trfBus(:,2:end),buslist(XfrmrKeep))); find(ismemberi(trfBus(:,3),buslist(XfrmrKeep)))];
	%find which transformers have delta connection and keep
% 	Delta_xfrmrs=find(sum(trf_conn,2)>0);
	Delta_xfrmrs=find(trfKeep);

	
	%find all transformers which should be kept. This is all Xfrmrs which have
	%delta connection and a CB behind them.
	XfrmrToKeep=XfrmrKeepInd(find(ismemberi(XfrmrKeepInd,Delta_xfrmrs)));
	
	%Remove transformers to keep from candidate to remove list
	candidateToRemove(find(ismember(candidateToRemove,XfrmrToKeep)))=[];
	
	%find regulator xfrmrs and keep
% 	nonRegInd=find(~ismemberi(1:length(circuit.transformer),Trf_Ind));
	%Remove from list of candidates
% 	candidateToRemove(find(~ismemberi(candidateToRemove,nonRegInd)))=[];
	candidateToRemove(find(ismemberi(candidateToRemove,Trf_Ind)))=[];

	%remove remaining transformers (not sub, reg, or delta with CB) from
	%circuit
	XfrmrToRmv=candidateToRemove;
	circuit.transformer(XfrmrToRmv)=[];
	
	rembus=[circuit.transformer{:}.buses]; rembus=unique(strtok(rembus(:),'\.'));
	rembus=find(ismemberi(buslist,rembus));
	criticalNumbers=[criticalNumbers; rembus];
	criticalNumbers=unique(criticalNumbers);
	criticalBuses=buslist(criticalNumbers);
	
	save(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Xfrmrs.mat'],'criticalNumbers','criticalBuses','circuit')
	clear trfDown trfBus trf_kV trfDownInd

else
	load(['c:\users\zactus\FeederReductionRepo\Feeders\' feeder '\' feeder '_Xfrmrs.mat'])
	end
end

t_=toc;
fprintf('%.2f sec\n',t_)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     TOPOGRAPHY NODES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this will now also become a critical node. This is done by looking at the
%grandparents of each critical node and seeing if there are common nodes
%between two critical nodes. The closest common node is the one that is
%kept
tic
fprintf('\nGetting topogrophical critical Nodes (connection): ')
nn=length(criticalNumbers);
% New_criticalNumbers1=[];
New_criticalNumbers2=[];


tic
if nn<2000
	for k=1:nn
		CB_parents=cell2mat(generation(criticalNumbers(k),4));
		c = generation(criticalNumbers(k+1:nn),4)';
		lens = sum(cellfun('length',c),1); innerCB_parents = ones(max(lens),numel(lens));
		innerCB_parents(bsxfun(@le,[1:max(lens)]',lens)) = vertcat(c{:});
		innerCB_parents(find(~ismemberi(innerCB_parents,CB_parents)))=1;
		VectorOfCommonParents=innerCB_parents;
		[~,num]=max(cell2mat(reshape(generation(VectorOfCommonParents,5),size(VectorOfCommonParents))));
		if ~isempty(num)
			New_criticalNumbers2=[New_criticalNumbers2,VectorOfCommonParents(sub2ind(size(VectorOfCommonParents), num, [1:length(num)]))];
		end
	end
else
	parfor k=1:nn
		CB_parents=cell2mat(generation(criticalNumbers(k),4));
		c = generation(criticalNumbers(k+1:nn),4)';
		lens = sum(cellfun('length',c),1); innerCB_parents = ones(max(lens),numel(lens));
		innerCB_parents(bsxfun(@le,[1:max(lens)]',lens)) = vertcat(c{:});
		innerCB_parents(find(~ismemberi(innerCB_parents,CB_parents)))=1;
		VectorOfCommonParents=innerCB_parents;
		[~,num]=max(cell2mat(reshape(generation(VectorOfCommonParents,5),size(VectorOfCommonParents))));
		if ~isempty(num)
			New_criticalNumbers2=[New_criticalNumbers2,VectorOfCommonParents(sub2ind(size(VectorOfCommonParents), num, [1:length(num)]))];
		end
	end
end
criticalNumbers=vertcat(criticalNumbers, unique(New_criticalNumbers2'));

%Add sourcebus to CN
criticalNumbers=[criticalNumbers; find(ismemberi(buslist,'sourcebus'))];

%Add all substation equipment to list of critical buses
criticalNumbers=[criticalNumbers; find(not(cellfun('isempty', strfind(lower(buslist), 'sub'))))];

criticalNumbers=unique(criticalNumbers); %get rid of repeat connections
criticalBuses=buslist(criticalNumbers);

t_=toc;
fprintf('time elapsed %f\n',t_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REDUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actual Reduction Part
fprintf('\n-------------------------------------------------------------\n')
fprintf('\n                     Reducing circuit          \n')
fprintf('\n-------------------------------------------------------------\n')

tic
fprintf('\nReducing: ')

YbusOrder_reduced=find(ismemberi(YbusOrderVect,criticalBuses));
Zbus_new=Zbus(YbusOrder_reduced,YbusOrder_reduced);
Ybus_reduced=inv(Zbus_new);
Ybus_reduced(find(abs(Ybus_reduced)<1E-7))=0;



%keep track of order
buslist_org=buslist; YbusOrderVect_org=YbusOrderVect; YbusPhaseVect_org=YbusPhaseVect; Ycomb_org=Ycomb;

buslist=criticalBuses; YbusOrderVect=YbusOrderVect(YbusOrder_reduced); YbusPhaseVect=YbusPhaseVect(YbusOrder_reduced); Ycomb=Ycomb(YbusOrder_reduced);

% W=Zbus_new\Zbus(YbusOrder_reduced,:);
W=Ybus_reduced*Zbus(YbusOrder_reduced,:);

W(find(abs(W)<1E-6))=0;


if Flag(1)
	S_LD_new=V2(YbusOrder_reduced,YbusOrder_reduced)*conj(W)*diag(1./diag(V2))*S_LD;
	fprintf('\nDiff between total load before and after reduction: %f kW\n',sum(sum(S_LD))-sum(sum(S_LD_new)))
end
if Flag(2)
	S_PV_new=V2(YbusOrder_reduced,YbusOrder_reduced)*conj(W)*diag(1./diag(V2))*S_PV;
	fprintf('\nDiff between pv before and after reduction: %f kW\n',sum(sum(S_PV))-sum(sum(S_PV_new)))
end

t_=toc;
fprintf('%.2f sec\n',t_)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REWRITE CIRCUIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-------------------------------------------------------------\n')
fprintf('\n                     Re-writing circuit        \n')
fprintf('\n-------------------------------------------------------------\n')

tic
fprintf('\nGetting Updated Topology: ')
%new topo
Node_number=[];
for ii=1:length(buslist)
	Ind=find(strcmpi(buslist(ii),YbusOrderVect))';
	Node_number(Ind)=ii;
end
topo=zeros(max(Node_number),4);
generation{1,1}=[];clear generation;  generation{1,1}=1; generation{1,4}=[];generation{1,5}=0;
parent=1;
topo(parent,1)=parent;
[topo,generation]=topology_detect_large(topo,generation,Ybus_reduced,parent,Node_number);
topo_view=topo;
topo_view(find(topo_view(:,1)==0)',:)=[];
c_new=0;

t_=toc;
fprintf('%.2f sec\n',t_)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   CLEAN-UP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names=fieldnames(circuit);
keepFields={'load','buslist','line','circuit','capcontrol','transformer','capacitor','basevoltages','regcontrol','pvsystem','reactor'};
names=names(find(~ismemberi(names,keepFields)));
for ii=1:length(names)
	circuit=rmfield(circuit,names{ii});
end
trf=circuit.transformer;
% for ii=1:length(trf)
% 	if([trf(ii).kv{1}]==696969)
% 		continue
% 	else
% 		Buses=trf(ii).Buses;
% 		for jj=1:length(Buses)
% 			trfBus(ii,jj)=lower(strtok(Buses(jj),'\.'));
% 			if trf(ii).phases==1
% 				trf_kV2(ii,jj)={cell2mat(trf(ii).kv(jj))*sqrt(3)};
% 			else
% 				trf_kV2(ii,jj)=trf(ii).kv(jj);
% 			end
% 		end
% 	end
% end
for ii=1:length(trf)
	if([trf(ii).kv{1}]==696969)
		continue
	else
		[bus,phase]=strtok(trf(ii).Buses,'.');
		trfBus(ii,1:length(bus))=bus;
		if trf(ii).phases==1
			trf_kV2(ii,1:length(trf(ii).kv))=cell2mat(trf(ii).kv)*sqrt(3);
		else
			trf_kV2(ii,1:length(trf(ii).kv))=cell2mat(trf(ii).kv);
		end
	end
end

%Get buses downstream of Vreg to allocate Vbase to load and PV
trfDown=trfBus(:,2:end); trfDown=trfDown(:); trfDown=trfDown(find(~cellfun(@isempty,trfDown)));
trf_kV=trf_kV2(:,2:end); trf_kV=trf_kV(:); %trf_kV=trf_kV(find(~cellfun(@isempty,trf_kV)));  trf_kV=cell2mat(trf_kV);
trf_kV(find(trf_kV==0))=[];
[~,trfDownInd]=ismemberi(trfDown,buslist);
%remove trf without kV specifications, likely just Vreg's
trfBus(cellfun(@isempty,trfBus))={'noBus'};

% [trfDownInd, IA]=unique(trfDownInd); trfDown=trfDown(IA);
% trf_kV=cell2mat(trf_kV(IA));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   LINES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
fprintf('\n Writing Lines: ')

%remove liens that are not enabled.
enabledBuses=find(~ismemberi(circuit.line(:).enabled,'false'));

Lines(:,1)=strtok(circuit.line(enabledBuses).bus1,'.');
Lines(:,2)=strtok(circuit.line(enabledBuses).bus2,'.');
Lines(:,3)=circuit.line(enabledBuses).length;
Lines(:,4)=circuit.line(enabledBuses).units;
[Lines(:,3),Lines(:,4)]=Convert_to_kft(Lines(:,4),Lines(:,3));

%find lines that are switches
Inds=[find(ismemberi(circuit.line(enabledBuses).switch,'true')) find(ismemberi(circuit.line(enabledBuses).switch,'yes'))];
Lines(Inds,3)={.001};

circuit.line=dssline;
count=0;
for ii=1:length(generation(:,1))
	Connected=cell2mat(generation(ii,2));
	Bus1=buslist(cell2mat(generation(ii,1)));
	Bus1Ind=find(ismemberi(YbusOrderVect,Bus1));
	
	%Make sure Phases are in correct order
	[~,I]=sort(YbusPhaseVect(Bus1Ind));
	Bus1Ind=Bus1Ind(I);
	
	for jj=1:length(Connected)
		Bus2=buslist(Connected(jj));
		Bus2Ind=find(ismemberi(YbusOrderVect,Bus2));
		[~,I]=sort(YbusPhaseVect(Bus2Ind));
		Bus2Ind=Bus2Ind(I);
		
		
		downSide=find(ismemberi(trfBus(:,1),Bus1));
		upSide=find(ismemberi(trfBus(:,2:end),Bus2));
		Match=find(ismemberi(downSide,upSide));
		%Make sure it is not VR
		if Match
			continue
		else
			count=count+1;
			
			circuit.line(count)=dssline;
			circuit.line(count).Name=[char(Bus1) '_' char(Bus2)];
			circuit.line(count).Units='kft';
			circuit.line(count).R1=[]; circuit.line(count).R0=[];
			circuit.line(count).X0=[]; circuit.line(count).X1=[];
			circuit.line(count).C0=[]; circuit.line(count).C1=[];
			
			Bus1IndMod=Bus1Ind(find(ismemberi(YbusPhaseVect(Bus1Ind),YbusPhaseVect(Bus2Ind))));
			circuit.line(count).bus1=[char(Bus1) '.' strjoin(arrayfun(@(x) num2str(x),YbusPhaseVect(Bus1IndMod),'UniformOutput',false),'.')];
			circuit.line(count).bus2=[char(Bus2) '.' strjoin(arrayfun(@(x) num2str(x),YbusPhaseVect(Bus2Ind),'UniformOutput',false),'.')];
			circuit.line(count).Phases=length(Bus2Ind);
			
			
			%calc length: Eseentially take bus 2 and back track to bus 1
			FinalLine=find(ismemberi(Lines(:,1),char(Bus1)));
			CurrentLine=find(ismemberi(Lines(:,2),char(Bus2)));
			LineInd=CurrentLine;
			Length=cell2mat(Lines(CurrentLine,3));
			while CurrentLine~=FinalLine
				CurrentLine=find(ismemberi(Lines(:,2),Lines(CurrentLine,1)));
				LineInd=[LineInd; CurrentLine];
				Length=Length+cell2mat(Lines(CurrentLine,3));
			end
			Lines(LineInd,:)=[];
			if ~isempty(Length)
				if length(Length)>1
					circuit.line(count).Length=mean(Length);%*0.3048;
					if mean(Length)~=Length(1)
						error('Logic is not correct here! Lines are different lengths!')
					end
				else
					circuit.line(count).Length=Length;%*0.3048;
				end
			else
				circuit.line(count).Length=.0001;
			end
			
			ImpMat=-inv(Ybus_reduced(Bus1IndMod,Bus2Ind))./(circuit.line(count).Length);
			ImpMatReal=real(ImpMat);
			ImpMatReal(abs(ImpMatReal)<1E-8)=0;
			ImpMatImag=imag(ImpMat);
			ImpMatImag(abs(ImpMatImag)<1E-8)=0;
			
			if length(ImpMat)==1
				circuit.line(count).Rmatrix=['(' num2str(ImpMatReal(1,1),6) ')'];
				circuit.line(count).Xmatrix=['(' num2str(ImpMatImag(1,1),6) ')'];
			elseif length(ImpMat)==2
				circuit.line(count).Rmatrix=['(' num2str(ImpMatReal(1,1),6) '|' num2str(ImpMatReal(2,1:2),6)  ')'];
				circuit.line(count).Xmatrix=['(' num2str(ImpMatImag(1,1),6) '|' num2str(ImpMatImag(2,1:2),6)  ')'];
			else
				circuit.line(count).Rmatrix=['(' num2str(ImpMatReal(1,1),6) '|' num2str(ImpMatReal(2,1:2),6) '|' num2str(ImpMatReal(3,1:3),6) ')'];
				circuit.line(count).Xmatrix=['(' num2str(ImpMatImag(1,1),6) '|' num2str(ImpMatImag(2,1:2),6) '|' num2str(ImpMatImag(3,1:3),6) ')'];
			end
			
		end
	end
end
t_=toc;
fprintf('%.2f sec\n',t_)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     CAPACITORS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%update capcontrol to match updated line names
if Flag(6)
	for ii=1:length(circuit.capcontrol)
		CapInd=find(ismemberi(circuit.capacitor(:).name,{circuit.capcontrol(ii).Capacitor}));
		bus1=strtok(circuit.capacitor(CapInd).bus1,'\.');
		bus2=strtok(circuit.capacitor(CapInd).bus2,'\.');
		
		if ~isempty(bus2)
			LineBus2=find(ismemberi(strtok(circuit.line(:).bus1,'.'),bus1));
			LineBus1=find(ismemberi(strtok(circuit.line(:).bus2,'.'),bus2));
			lineInd=LineBus2(find(ismemberi(LineBus2,LineBus1)));
		else
			nextBus=generation{find(ismemberi(buslist,bus1)),2};
			prevBus=generation{find(ismemberi(buslist,bus1)),4};
			
			%if line exists to bus after cap bus, use that line, else use
			%line to previous bus.
			if isempty(nextBus)
				LineBus2=find(ismemberi(strtok(circuit.line(:).bus2,'.'),bus1));
				LineBus1=find(ismemberi(strtok(circuit.line(:).bus1,'.'),buslist(prevBus(1))));
			else
				LineBus1=find(ismemberi(strtok(circuit.line(:).bus1,'.'),bus1));
				LineBus2=find(ismemberi(strtok(circuit.line(:).bus2,'.'),buslist(nextBus(1))));
			end
			lineInd=LineBus2(find(ismemberi(LineBus2,LineBus1)));
		end
		
		circuit.capcontrol(ii).Element=['line.' circuit.line(lineInd).Name];
	end
end

%add shunt cap
%find rows with non-zero shunt capacitanct
ShuntRows=find(sum(Ybus_reduced,2)>1E-6);

% 	if any(ismemberi(trfIndYbus,keepRows(ii)))
% 		%get Yprim and subtract from row
% 		[trfInd,~]=ind2sub(size(trfBus),find(ismemberi(trfBus,YbusOrderVect(keepRows(ii)))));
% 		dssCircuit.SetActiveElement(['transformer.' circuit.transformer(unique(trfInd)).name]);
% 		Yprim=dssCircuit.ActiveElement.Yprim;
%         Yprim0=Yprim(1:2:end)+1i*Yprim(2:2:end);
%         NewYprimLength=sqrt(length(Yprim0));
%         Yprim=reshape(Yprim0,[NewYprimLength,NewYprimLength]);
%
% 	end

%find rows that correspond to a transformer in the system
trfIndYbus=[find(ismemberi(Ycomb,regexprep([circuit.transformer{:}.buses],'\.0',''))); find(ismemberi(YbusOrderVect,regexprep([circuit.transformer{:}.buses],'\.0','')))];
% ShuntRows=ShuntRows(find(~ismember(ShuntRows,trfIndYbus)));

%remove rows that correspond to a capacitor in the system
busCap=[circuit.capacitor(:).bus1; circuit.capacitor(:).bus2];
busCap=busCap(find(~cellfun(@isempty,busCap)));
capIndYbus=[find(ismemberi(Ycomb,busCap)); find(ismemberi(YbusOrderVect,busCap))];
ShuntRows=ShuntRows(find(~ismember(ShuntRows,capIndYbus)));

%remove sourcebus
Source=find(ismemberi(YbusOrderVect,'sourcebus'));
keepRows=ShuntRows(find(~ismember(ShuntRows,Source)));
%  keepRows=[];

% %Make additional capacitors here to account for removed impedence from DT
for ii=1:length(keepRows)
	if any(ismemberi(trfIndYbus,keepRows(ii)))
		r_jx=1/(sum(Ybus_reduced(keepRows(ii),:),2)-sum(YofXfrmrs(find(ismemberi(Ycomb_org,Ycomb(keepRows(ii)))),:)));
	else
		r_jx=1/sum(Ybus_reduced(keepRows(ii),:),2);
	end
	
	if real(1/r_jx)<1E-6 && imag(1/r_jx)<1E-6
		continue
	elseif real(1/r_jx)<1E-6
		r_jx=0+1i*imag(r_jx);
	elseif imag(1/r_jx)<1E-6
		r_jx=real(r_jx)+0;
	end
	
	if ~Flag(7) && ii==1
		circuit.reactor(1)=dssreactor;
	else
		circuit.reactor(end+1)=dssreactor;
	end
	
	
	circuit.reactor(end).name=['addedShunt_' char(regexprep(Ycomb{keepRows(ii)},'\.','_'))];
	circuit.reactor(end).phases=1;
	circuit.reactor(end).bus1=char(Ycomb(keepRows(ii)));
	
	circuit.reactor(end).R=real(r_jx);
	circuit.reactor(end).X=imag(r_jx);
	circuit.reactor(end).kVAr=0;
	
	busNum=find(ismemberi(buslist,YbusOrderVect(keepRows(ii))));
	MatchInd=find(ismemberi(trfDownInd,[cell2mat(generation(busNum,4)); busNum]));
	
	circuit.reactor(end).kv=trf_kV(MatchInd(end))/sqrt(3);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     LOADS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Flag(1)
	fprintf('Re-writing loads: ')
	tic
	
	count=0;
	circuit.load=dssload;
	for loadModel=1:size(S_LD_new,2)
		WriteLoads=find(abs(S_LD_new(:,loadModel))>0.01);
		for ii=1:length(WriteLoads)
			count=count+1;
			circuit.load(count)=dssload;
			circuit.load(count).Name=['Load ' char(Ycomb(WriteLoads(ii))) '_model_' num2str(loadModel)];
			circuit.load(count).phases=1;
			
			busNum=find(ismemberi(buslist,YbusOrderVect(WriteLoads(ii))));
			MatchInd=find(ismemberi(trfDownInd,[cell2mat(generation(busNum,4)); busNum]));
			circuit.load(count).kV= trf_kV(MatchInd(end))/sqrt(3);
			circuit.load(count).bus1=Ycomb(WriteLoads(ii));
			circuit.load(count).Kw=real(S_LD_new(WriteLoads(ii),loadModel));
			circuit.load(count).KVAr=imag(S_LD_new(WriteLoads(ii),loadModel));
			circuit.load(count).model=loadModel;
			circuit.load(count).kVA=[];
			circuit.load(count).pf=[];
		end
	end
	
	if strcmpi(feeder,'IEEE8500')
		for ii=1:length(circuit.load)
			circuit.load(ii).vminpu=.88;
		end
	end
	
	t_=toc;
	fprintf('%.2f sec\n',t_)
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PV    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Flag(2)
	fprintf('Re-writing PV systems: ')
	tic
	
	count=0;
	circuit.pvsystem=dsspvsystem;
	WritePVs=find(abs(S_PV_new(:))>0.01);
	for ii=1:length(WritePVs)
		count=count+1;
		circuit.pvsystem(count)=dsspvsystem;
		circuit.pvsystem(count).Name=['PV_' char(Ycomb(WritePVs(ii)))];
		circuit.pvsystem(count).phases=1;
		circuit.pvsystem(count).irradiance=1;
		% 		circuit.pvsystem(count).Temperature=25;
		
		busNum=find(ismemberi(buslist,YbusOrderVect(WritePVs(ii))));
		MatchInd=find(ismemberi(trfDownInd,[cell2mat(generation(busNum,4)); busNum]));
		circuit.pvsystem(count).kV=trf_kV(MatchInd(end))/sqrt(3);
		trf_kV(MatchInd(end))/sqrt(3);
		
		if strcmp(circuit.pvsystem(count).bus1,'03551328A.1')
			stop=1;
		end
		
		circuit.pvsystem(count).bus1=Ycomb(WritePVs(ii));
		circuit.pvsystem(count).pmpp=sqrt(real(S_PV_new(WritePVs(ii)))^2+imag(S_PV_new(WritePVs(ii)))^2);
		circuit.pvsystem(count).kVA=sqrt(real(S_PV_new(WritePVs(ii)))^2+imag(S_PV_new(WritePVs(ii)))^2);
		% 		circuit.pvsystem(count).pf=real(S_PV_new(WritePVs(ii)))/circuit.pvsystem(count).kVA;
		%Fix PF sign
		if sign(real(S_PV_new(ii)))*sign(imag(S_PV_new(ii)))<0
			% 			circuit.pvsystem(count).pf=-circuit.pvsystem(count).pf;
		end
		circuit.pvsystem(count).cutout=0;
		circuit.pvsystem(count).cutin=0;
		
	end
	
	t_=toc;
	fprintf('%.2f sec\n',t_)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Update circuit info
circuit.circuit.Name=[circuit.circuit.Name '_Reduced'];
tmp=find(ismemberi(circuit.buslist.id,buslist));
circuit.buslist.id=circuit.buslist.id(tmp);
circuit.buslist.coord=circuit.buslist.coord(tmp,:);

%store Ybus
circuit.Ybus=Ybus_reduced;
circuit.YbusOrder=Ycomb;

t_reduct=toc(Full);
if Flag(2)
	circuitBase=rmfield(circuit,'pvsystem');
else
	circuitBase=circuit;
end
if Flag(1);
	circuitBase=rmfield(circuitBase,'load');
end
if Flag(5)
	circuitBase=rmfield(circuitBase,'regcontrol');
end

p = WriteDSS(circuitBase,'Test',0,'c:\users\zactus\FeederReductionRepo\Feeders\IEEE8500\tmp');
% p = WriteDSS(circuit,'Test',0,'c:\users\zactus\FeederReductionRepo\Feeders\IEEE8500\tmp');
fprintf('total reduction: %.2f sec\n',t_reduct)

[YbusOrderVect2, YbusPhaseVect2, Ycomb2, Ybus_regenerated, ~,~]=getYbus(circuitBase);

%plotting
redd=tic;
[ data ] = dssSimulation( circuit,[],1,[],[],0);
redTime=toc(redd);

for ii=1:length(data.nodeName)
	Keep(ii)=find(ismemberi(data1.nodeName,data.nodeName(ii)));
end

[diffV, ind]=max(abs(data1.Voltage(Keep)-data.Voltage));
diffV
Ycomb(ind)

figure;plot(data.Dist*.3048,data.Voltage,'r*',data1.Dist(Keep),data1.Voltage(Keep),'b*')
legend('Reduced','Original')
title('Matching Buses','fontsize',14)
xlabel('Distance to substation [km]','fontsize',14)
ylabel('Bus Voltage [V pu]','fontsize',14)
figure;plot(data.Dist*.3048,data.Voltage,'r*',data1.Dist,data1.Voltage,'b*')
legend('Reduced','Original')
title('All Buses','fontsize',14)
xlabel('Distance to substation [km]','fontsize',14)
ylabel('Bus Voltage [V pu]','fontsize',14)
figure;plot(1:length(Keep),data.Voltage,'r*',1:length(Keep),data1.Voltage(Keep),'b*')
legend('Reduced','Original')
title('Matching Buses','fontsize',14)
xlabel('Bus number [-]','fontsize',14)
ylabel('Bus Voltage [V pu]','fontsize',14)

keySet = lower(Ycomb2);
valueSet = lower(Ycomb);
mapObj = containers.Map(keySet,1:length(keySet));
Matching_orig_order=cell2mat(values(mapObj,valueSet));
Ybus_regenerated=Ybus_regenerated(Matching_orig_order,Matching_orig_order);
Ydiff=Ybus_regenerated-Ybus_reduced;
Ydiff(find(abs(Ydiff)<1E-3))=0;
[maxYbusDiff,ind] = max(abs(Ydiff(:)));
if abs(maxYbusDiff)>0
	maxYbusDiff
	[m,n] = ind2sub(size(Ydiff),ind);
	YbusOrderVect([m n])
end



% function [Output]=Convert_to_kft(units, Input)
%
% if strcmp(units,'kft')
% 	Output=Input;
% elseif strcmp(units,'mi')
% 	Output=Input*5.280;
% elseif strcmp(units,'km')
% 	Output=Input*3.28084;
% elseif strcmp(units,'m')
% 	Output=Input*.00328084;
% elseif strcmp(units,'ft')
% 	Output=Input/1000;
% elseif strcmp(units,'in')
% 	Output=Input/12/1000;
% elseif strcmp(units,'cm')
% 	Output=Input/2.54/12/1000;
% elseif strcmp(units,'none')
% 	Output=Input;
% end
% end






