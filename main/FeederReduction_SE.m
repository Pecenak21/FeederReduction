function [circuit,outputdss] = FeederReduction_SE(Critical_buses,circuit, Measurements)

%Created by Zachary K. Pecenak on 11/2/2016

%The purpose of this function is to reduce an OpenDSS feeder to a few
%selected node for state estimation. A user might want to reduce a feeder for several reasons,
%but mainly to reduce the system complexity for increased computational
%efficiency. This function works by taking in the openDSS feeder model and
%reducing the impedences, loads, and PV to a few specified nodes in terms of measurements. Note that
%nodes connecting two or more critical nodes are inherently also critical
%nodes.

%The critical nodes are the ones that the user wants to keep in the reduced
%circuit. Thdey are input as a cell of the node names in the circuit file
%or the node number

%The circuit is the c file for the specific feeder setup, and can be found
%in the saved feeder setup file.
%weights are stored in circuit.weights

%The measurments are given as a matrix and have the followign form
%accordign to Vahid R. Disfani's DSSE code which is the code this code was
%written for.
% M1: Pinj measurements     {1, full node name, 0, value, standard deviation}
% M2: Qinj measurements     {2, full node name, 0, value, standard deviation}
% M3: Pline measurements    {3, sending node name, receiving node name, value, standard deviation}
% M4: Qline measurements    {4, sending node name, receiving node name, value, standard deviation}
% M5: Voltage measurements  {5, full node name, 0, value, standard deviation}

%Example input
% Critical_buses={'03551325','03553322','035560','03555704','0355948','03552613','03552637','03551382','03554401','03552641A'};
% [c] = FeederReduction(Critical_buses,c);





%% TODO
%1) Fix ES
%2) Fix capacitive load
%3) Test All changes made
%4) Figure out y-delta connections for transformers!!!!
%5) Make transformers a line

%Check both inputs are met
if nargin<2
	error('This requires two inputs, the desired nodes to keep and the feeder circuit')
end
keepTrans=0;

Batt=[];
PV=[];
LOAD=[];
Weights=[];
InjP=[];
InjQ=[];
MeasV=[];
LineP=[];
LineQ=[];

pvFlag=0;
battFlag=0;
loadFlag=0;
MeasFlag=0;

if isfield(circuit,'pvsystem')
	pvFlag=1;
end
if isfield(circuit,'load')
	loadFlag=1;
end
if isfield(circuit,'storage')
	battFlag=1;
end

if nargin>2
	MeasFlag=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial cleanup of data

if MeasFlag
	% Fix Z to make sure that all nodes are represented and split into seperate
	% vectors for ease of programming
	Pinj=Measurements(find([Measurements{:,1}]==1),:);
	Qinj=Measurements(find([Measurements{:,1}]==2),:);
	Pline=Measurements(find([Measurements{:,1}]==3),:);
	Qline=Measurements(find([Measurements{:,1}]==4),:);
	Vmeas=Measurements(find([Measurements{:,1}]==5),:);
	
	clear Measurements
end
%make sure nodes are column
if ~iscolumn(Critical_buses)
	Critical_buses=Critical_buses';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get Ybus from OpenDSS for feeder.
fprintf('\nGetting YBUS and buslist from OpenDSS: ')
tic

if ~exist(['c:\users\zactus\gridIntegration\NewRes\' circuit.circuit.Name '_Ybus.mat'])
	%remove load and PV from Ybus
	if isfield(circuit,'pvsystem')
		circuit_woPV=rmfield(circuit,'pvsystem');
	else
		circuit_woPV=circuit;
	end
	if isfield(circuit,'load');
		circuit_woPV=rmfield(circuit_woPV,'load');
	end
	if isfield(circuit,'storage');
		circuit_woPV=rmfield(circuit_woPV,'storage');
	end
	%load the circuit and generate the YBUS
	p = WriteDSS(circuit_woPV,[],0,[]); o = actxserver('OpendssEngine.dss');
	dssText = o.Text; dssText.Command = 'Clear'; cDir = pwd;
	dssText.Command = ['Compile "' p '"']; dssCircuit = o.ActiveCircuit;
	Ybus=dssCircuit.SystemY;
	
	%Convert the Ybus to a matrix
	ineven=2:2:length(Ybus); inodd=1:2:length(Ybus);
	Ybus=Ybus(inodd)+1i*Ybus(ineven); Ybus=reshape(Ybus,sqrt(length(Ybus)),sqrt(length(Ybus)));
	Ybus=sparse(Ybus);
	%get buslist in order of Ybus and rearrange
	busnames=regexp(dssCircuit.YNodeOrder,'\.','split');
	YbusOrderVect=[busnames{:}]'; YbusOrderVect(find(cellfun('length',YbusOrderVect)==1))=[];
	YbusPhaseVect=[busnames{:}]'; YbusPhaseVect(find(cellfun('length',YbusPhaseVect)>1))=[]; YbusPhaseVect=str2double(YbusPhaseVect);
	Ycomb=strcat(YbusOrderVect,'.', num2str(YbusPhaseVect));
	% %Here we generate the list of node numbers in the order of the Ybus
	buslist=dssCircuit.AllBusNames;
	Origbuslist=buslist;
	for ii=1:length(buslist)
		Ind=find(strcmpi(buslist(ii),YbusOrderVect))';
		Node_number(Ind)=ii;
	end
	
	clear inodd ineven
	save(['c:\users\zactus\gridIntegration\NewRes\' circuit.circuit.Name '_Ybus.mat'],'Ycomb','Node_number','Origbuslist','buslist','YbusPhaseVect','YbusOrderVect','busnames','Ybus')
	delete(o)
else
	load(['c:\users\zactus\gridIntegration\NewRes\' circuit.circuit.Name '_Ybus.mat'])
end

t_=toc;
fprintf('time elapsed %f',t_)

%Check to see that the critical nodes are in the circuit

% Critical_buses(find(~ismember(lower(Critical_buses),lower(buslist))))=[];
% if ~any(ismember(lower(Critical_buses),lower(buslist)))
%  	error('The following nodes arent in the circuit: \n%s', Critical_buses{~ismember(Critical_buses,buslist)})
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get topogrophy %% written by Vahid R. Disfani
%Uncommented
fprintf('\nGetting topogrophy: ')
tic

if ~exist(['c:\users\zactus\gridIntegration\NewRes\' circuit.circuit.Name '_OrigTop.mat'])
	
	topo=zeros(max(Node_number),4);
	generation{1,1}=[];clear generation;  generation{1,1}=1; generation{1,4}=[];generation{1,5}=0;
	parent=1;
	topo(parent,1)=parent;
	[topo,generation]=topology_detect_large(topo,generation,Ybus,parent,Node_number);
	topo_view=topo;
	topo_view(find(topo_view(:,1)==0)',:)=[];
	c_new=0;
	
	save(['c:\users\zactus\gridIntegration\NewRes\' circuit.circuit.Name '_OrigTop.mat'],'generation','topo','topo_view')
	
else
	load(['c:\users\zactus\gridIntegration\NewRes\' circuit.circuit.Name '_OrigTop.mat'])
	c_new=0;
end
t_=toc;
fprintf('time elapsed %f',t_)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get line lengths and buses
fprintf('\nAdding lines: ')
tic

if ~exist(['c:\users\zactus\FeederReductionRepo\' circuit.circuit.Name '_Lines.mat'])
	%account for units, convert all to km
	I_cap=zeros(length(YbusOrderVect),1);
	jw=2*pi*60;
	
	trf_bus=strtok(reshape([circuit.transformer{:}.buses],2,[]),'\.'); trf_bus(1,:)=[];
	[~,trf_bus_ind]=ismember(lower(trf_bus),lower(buslist));
	[trf_bus_ind, IA]=unique(trf_bus_ind); trf_bus=trf_bus(IA);
	trf_kV=cell2mat(reshape([circuit.transformer{:}.Kv],2,[])); trf_kV(1,:)=[];
	trf_phase=[circuit.transformer{:}.Phases];
	trf_phase(find(trf_phase==3))=4; trf_phase(find(trf_phase==1))=3; trf_phase(find(trf_phase==4))=1;
	trf_kV=trf_kV.*sqrt(trf_phase);  trf_kV=trf_kV(IA);
	
	C_MAT=zeros(length(YbusOrderVect));
	
	for ii=1:length(circuit.line)
		norm_flag=0;
		LnCd_ind=find(ismember(lower({circuit.linecode{:}.Name}),lower(circuit.line(ii).LineCode)));
		
		%Convert all lengths to kft
		units=lower(circuit.line(ii).Units);

		if isempty(units)
			units=lower(circuit.linecode(LnCd_ind).Units);
		end
		if strcmp(units,'none')
			norm_flag=1; %Stupid OpenDSS. This accounts for the fact that impedence is not per unit length
		end
		
		[Lines{ii,1}]=Convert_to_kft(units, circuit.line(ii).Length);
		
		Lines{ii,2}=lower(regexprep(circuit.line(ii).bus1,'\..*',''));
		Lines{ii,3}=lower(regexprep(circuit.line(ii).bus2,'\..*',''));
		
		if ~isempty(circuit.line(ii).LineCode)
			units=lower(circuit.linecode(LnCd_ind).Units);
			if ~isempty(circuit.linecode(LnCd_ind).Cmatrix)
				if ~ischar(circuit.linecode(LnCd_ind).Cmatrix)
					Cmat=Convert_to_kft(units,circuit.linecode(LnCd_ind).Cmatrix); %nano A/V/Length
				else
					c_tmp=regexprep(circuit.linecode(LnCd_ind).Cmatrix,'[','');
					c_tmp=regexprep(c_tmp,']','');
					c_tmp=regexp(c_tmp,'\|','split');
					for iiii=1:length(c_tmp)
						Cmat(iiii,1:iiii)=Convert_to_kft(units,str2num(cell2mat(c_tmp(iiii))));
					end
					for iiii=1:size(Cmat,1)
						for jjjj=iiii+1:size(Cmat,2)
							Cmat(iiii,jjjj)=Cmat(jjjj,iiii);
						end
					end
				end
			elseif ~isempty(circuit.linecode(LnCd_ind).C0)
				Cm=(1/3)*(circuit.linecode(LnCd_ind).C0-circuit.linecode(LnCd_ind).C1);
				Cs=circuit.linecode(LnCd_ind).C1+Cm;
				Cmat=Convert_to_kft(units,diag(Cs.*ones(circuit.line(ii).Phases,1))+Cm*ones(circuit.line(ii).Phases)-diag(Cm.*ones(circuit.line(ii).Phases,1)));
			else
				Cmat=zeros(circuit.line(ii).Phases);
			end
		else
			if ~isempty(circuit.line(ii).Cmatrix)
				if ~ischar(circuit.linecode(LnCd_ind).Cmatrix)
					Cmat=Convert_to_kft(units,circuit.line(ii).Cmatrix); %nano V/C/Length
				else
					c_tmp=regexprep(circuit.line(ii).Cmatrix,'[','');
					c_tmp=regexprep(c_tmp,']','');
					c_tmp=regexp(c_tmp,'\|','split');
					for iiii=1:length(c_tmp)
						Cmat(iiii,1:iiii)=Convert_to_kft(units,str2num(cell2mat(c_tmp(iiii))));
					end
					for iiii=1:size(Cmat,1)
						for jjjj=iiii+1:size(Cmat,2)
							Cmat(iiii,jjjj)=Cmat(jjjj,iiii);
						end
					end
				end
			elseif ~isempty(circuit.line(ii).C0)
				Cm=(1/3)*(circuit.line(ii).C0-circuit.line(ii).C1);
				Cs=circuit.line(ii).C1+Cm;
				Cmat=Convert_to_kft(units,diag(Cs.*ones(circuit.line(ii).Phases,1))+Cm*ones(circuit.line(ii).Phases)-diag(Cm.*ones(circuit.line(ii).Phases,1)));
			else
				Cmat=zeros(circuit.line(ii).Phases);
			end
		end
		
		[~,phase1]=strtok(circuit.line(ii).bus1,'\.');
		Phases1=str2num(regexprep(phase1,'\.','')');
		Inds_all1=find(ismember(lower(YbusOrderVect),lower(Lines{ii,2})));
		if isempty(phase1)
			Inds1=Inds_all1;
		else
			Inds1=Inds_all1(find(ismember(YbusPhaseVect(Inds_all1),Phases1)));
		end
		
		[~,phase2]=strtok(circuit.line(ii).bus2,'\.');
		Phases2=str2num(regexprep(phase2,'\.','')');
		Inds_all2=find(ismember(lower(YbusOrderVect),lower(Lines{ii,3})));
		if isempty(phase2)
			Inds2=Inds_all2;
		else
			Inds2=Inds_all2(find(ismember(YbusPhaseVect(Inds_all2),Phases2)));
		end
		
		%Special case for switch=true
		if strcmpi(circuit.line(ii).switch,'true') | strcmpi(circuit.line(ii).switch,'yes')
			Lines{ii,1}=.001;
			Cmat=zeros(circuit.line(ii).Phases);
		end
		%Multiply the cap matrix by the voltage and divide between phases:
		%I=VCjw
		
		
		%Make sure kV is correct
		BusWithLoad=find(ismember(buslist,lower(YbusOrderVect(Inds1))));
		MatchInd=find(ismember(trf_bus_ind,[cell2mat(generation(BusWithLoad,4)); BusWithLoad]));
			
		if norm_flag~=1;
			I_cap(Inds1)=I_cap(Inds1)+trf_kV(MatchInd(end))*1E-9*Cmat.*Lines{ii,1}*ones(length(Inds1),1).*jw./length(Inds1)/2;
			I_cap(Inds2)=I_cap(Inds2)+trf_kV(MatchInd(end))*1E-9*Cmat.*Lines{ii,1}*ones(length(Inds2),1).*jw./length(Inds2)/2;

		else
			I_cap(Inds1)=I_cap(Inds1)+trf_kV(MatchInd(end))*1E-9*Cmat*ones(length(Inds1),1).*jw./length(Inds1)/2;
			I_cap(Inds2)=I_cap(Inds2)+trf_kV(MatchInd(end))*1E-9*Cmat*ones(length(Inds2),1).*jw./length(Inds2)/2;
		end
		
		C_MAT(Inds1,Inds2)=Cmat;
	end
	
	save(['c:\users\zactus\FeederReductionRepo\' circuit.circuit.Name '_Lines.mat'],'C_MAT','I_cap','Lines')
	clear trf_bus_ind Cm Cmat
else
	load(['c:\users\zactus\FeederReductionRepo\' circuit.circuit.Name '_Lines.mat'])
end
t_=toc;
fprintf('time elapsed %f',t_)

%% Re-organize Measurements

if MeasFlag
	fprintf('\nOrganizing Measurements: ')
	tic
	%Now that we have node information and lien information, lets organize the data
	InjP=zeros(length(YbusOrderVect),2);
	InjQ=zeros(length(YbusOrderVect),2);
	MeasV=zeros(length(YbusOrderVect),2);
	
	LineP=Pline(:,2:5);
	LineQ=Qline(:,2:5);
	
	
	%Put data into vector arrays
	%Pinj
	for ii=1:length(Pinj)
		InjP(find(ismember(lower(Ycomb),lower({Pinj{ii,2}}))),:)=[Pinj{ii,4:5}];
	end
	%Qinj
	for ii=1:length(Qinj)
		InjQ(find(ismember(lower(Ycomb),lower({Qinj{ii,2}}))),:)=[Qinj{ii,4:5}];
	end
	%Vmeas
	for ii=1:length(Vmeas)
		MeasV(find(ismember(lower(Ycomb),lower({Vmeas{ii,2}}))),:)=[Vmeas{ii,4:5}];
	end
	
	t_=toc;
	fprintf('time elapsed %f',t_)
end

%% get pv data and match with nodes
%Format the various PV things
if isfield(circuit,'pvsystem')
	fprintf('\nGetting PV of circuit: ')
	if ~exist(['c:\users\zactus\FeederReductionRepo\' circuit.circuit.Name '_PV.mat'])
		
		%Allocate space
		pv=circuit.pvsystem;
		PV=zeros(length(YbusOrderVect),1);
		PVbusMap=zeros(length(pv),length(buslist));
		Weights=zeros(length(YbusOrderVect),length(buslist));
		PvSize=zeros(length(pv),1)';
		
		
		
		for j=1:length(pv)
		
			%Break up bus name to get bus and phases
			PVname=regexp(char(pv.bus1(j)),'\.','split','once');
			if length(PVname)>1
				PvPhase{j}=PVname{2};
			else
				PvPhase{j}='1.2.3';
			end
			
			%Get the bus and nodes of the PV
			PVbusInd(j)=find(ismember(lower(buslist),lower(PVname{1})));
			PVNodeInd=find(ismember(lower(YbusOrderVect),lower(PVname{1})));
			
			%Adjust Nodes to get accurate phases
			PvPhasesVect=regexp(PvPhase{j},'\.','split'); PvPhasesVect=str2num(cell2mat(PvPhasesVect(:)));
			
			%Assign correct phases to PV and update
			MatchInd=find(ismember(YbusPhaseVect(PVNodeInd),YbusPhaseVect(PvPhasesVect)));
			PV(PVNodeInd(MatchInd))=PV(PVNodeInd(MatchInd))+pv(j).kVa/length(MatchInd);
			
			%get weights %basically maps PV to node phases and buses, but gets
			%updataed later
			Weights(PVNodeInd(MatchInd),PVbusInd(j))=1/length(PVNodeInd(MatchInd));
			
			%Map PV to appropriate buses for later conversion
			PVbusMap(j,PVbusInd(j))=1;
			
			% Store vector of PV size for later conversion
			PvSize(j)=double(pv(j).kVa);
		end
		clearvars PVNodeInd PvPhasesVect PVbusInd PvPhase
		save(['c:\users\zactus\FeederReductionRepo\' circuit.circuit.Name '_PV.mat'],'PV','Weights','PVbusMap','PvSize')
	else
		load(['c:\users\zactus\FeederReductionRepo\' circuit.circuit.Name '_PV.mat'])
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get load data and match with nodes
%All of this is basically to just make the load 3 phase
%There is undoubtably a more efficient way
if isfield(circuit,'load')
	fprintf('\nGetting Load of circuit: ')
	if ~exist(['c:\users\zactus\FeederReductionRepo\' circuit.circuit.Name '_Load.mat'])
		
		ld=circuit.load;
		LOAD=zeros(length(YbusOrderVect),1);
		
		for j=1:length(ld)
			if isempty(regexp(ld(j).bus1,'\.','match'))
				LOAD(find(strcmpi(YbusOrderVect,ld(j).bus1)==1))=ld(j).Kw/3+i*ld(j).Kvar/3;
			elseif length(regexp(ld(j).bus1,'\.','match'))>1
				name=regexp(ld(j).bus1,'\.','split');
				for ii=2:length(name)
					LOAD(find(strcmpi(Ycomb,[name{1} '.' name{ii}])==1))=(ld(j).Kw/length(2:length(name))+i*ld(j).Kvar/length(2:length(name)));
				end
			else
				LOAD(find(strcmpi(Ycomb,ld(j).bus1)==1))=ld(j).Kw+i*ld(j).Kvar;
			end
		end
		save(['c:\users\zactus\FeederReductionRepo\' circuit.circuit.Name '_Load.mat'],'LOAD','Ycomb')
	else
		load(['c:\users\zactus\FeederReductionRepo\' circuit.circuit.Name '_Load.mat'])
	end
end

t_=toc;
fprintf('time elapsed %f',t_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ld_a_beg=sum(LOAD(find(YbusPhaseVect==1)));
Ld_b_beg=sum(LOAD(find(YbusPhaseVect==2)));
Ld_c_beg=sum(LOAD(find(YbusPhaseVect==3)));
%% Get Storage and such
%Things to consider:
%KV all the same?
%Aggregate KVA
%Aggregate kvar
%aggregate kWhrated
if isfield(circuit,'storage')
	tic
	fprintf('\nAdding storage devices: ')
	
	
	%Batt
	%1: kwhrated
	%2: kwhstored
	%3: kwrated
	%4: kva
	
	
	Batt=ones(length(YbusOrderVect),4);
	
	
	batt=circuit.storage;
	for ii=1:length(batt)
		
		%get bus and phase of battery
		[bus1,phase]=strtok(batt(ii).bus1,'\.');
		
		if length(phase)>1
			%split phase completely
			phase=regexprep(phase,'\.',''); phase=str2num(phase');
		elseif isempty(phase)
			phase=[1 2 3]';
		else
			phase=num2str(phase);
		end
		
		%find Indices
		BusInd=find(ismember(YbusOrderVect,bus1));
		BusInd=find(ismember(YbusPhaseVect(BusInd),phase));
		
		% 		Batt(BusInd,1)=batt.kwhrateded/length(BusInd);
		% 		Batt(BusInd,2)=batt.kwhstored/length(BusInd);
		% 		Batt(BusInd,3)=batt.kwrated/length(BusInd);
		% 		Batt(BusInd,4)=batt.kva/length(BusInd);
		
	end
	if isfield(circuit,'storagecontroller')
		battCont=circuit.storagecontroller;
		for ii=1:length(battCont)
			
			
			
		end
	end
	t_=toc;
	fprintf('time elapsed %f',t_)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% match critical nodes with node numbers
for j=1:length(Critical_buses)
	Critical_numbers(j)=find(strcmpi(buslist,lower(Critical_buses(j))));
end
Orig_nodes=Critical_buses;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add capacitors + transformers to Critical Nodes
%This section is to keep the capacitors and the transformers in the grid by
%adding the nodes that they are connected to.
bus1=[];
%Capacitor
if isfield(circuit,'capacitor')
	tic
	fprintf('\nAdding Capacitor Nodes: ')
	
	cap=circuit.capacitor;
	for ii=1:length(cap)
		bus1{ii}=cap(ii).bus1;
	end
	if ~iscolumn(bus1)
		bus1=bus1';
	end
	bus1=regexp(bus1,'\.','split');
	bus1=[bus1{:}]'; bus1(find(cellfun('length',bus1)==1))=[];
	New_Node_Num=find(ismember(lower(buslist),lower(bus1)));
	Critical_buses=[Critical_buses; lower(buslist(find(ismember(lower(buslist),lower(bus1)))))];
	Critical_numbers=[Critical_numbers'; New_Node_Num];
end
clear bus1

if isfield(circuit,'capcontrol')
	capcon=circuit.capcontrol;
	for ii=1:length(capcon)
		buses=regexp(capcon(ii).Element,'\.','split');
		if strcmp(buses{1},'line')
			LineNo=find(ismember({circuit.line{:}.Name},buses{2}));
			bus1{ii}=strtok(circuit.line(LineNo).bus1,'.');
			bus2{ii}=strtok(circuit.line(LineNo).bus2,'.');
		end
		%This was used to organize element by splitting name, stupid
		%Naming convention used in Kleissl models (i.e. line.bus1_bus2)
		% 		buses=regexp(buses{2},'\_','split');
		% 			bus1{ii}=buses{1};
		% 			bus2{ii}=buses{2};
	end
	% 	bus1=regexp(bus1,'\.','split');
	% 	bus1=[bus1{:}]'; bus1(find(cellfun('length',bus1)==1))=[];
	capConBuses1=bus1';
	% 	bus2=regexp(bus2,'\.','split');
	% 	bus2=[bus2{:}]'; bus2(find(cellfun('length',bus2)==1))=[];
	capConBuses2=bus2';
	New_Node_Num=reshape([find(ismember(lower(buslist),lower(bus1))), find(ismember(lower(buslist),lower(bus2)))],[],1);
	Critical_buses=[Critical_buses; lower(buslist(find(ismember(lower(buslist),lower(bus1))))); lower(buslist(find(ismember(lower(buslist),lower(bus2)))))];
	Critical_numbers=[Critical_numbers; New_Node_Num];
end
clear bus1 bus2
t_=toc;
fprintf('time elapsed %f',t_)
%% Remove unnecessary Xfrms

%Set all VR to CB (for now)
if isfield(circuit,'transformer')
	
	if isfield(circuit,'regcontrol')
		rgc=circuit.regcontrol;
		for ii=1:length(rgc)
			trfName=rgc(ii).transformer;
			buses=circuit.transformer(find(ismember({circuit.transformer{:}.Name},trfName))).buses;
			busNames=strtok(buses,'\.');
			Critical_numbers=[Critical_numbers; find(ismember(lower(buslist),lower(busNames)))];
		end
	end
	Critical_numbers=unique(Critical_numbers);
	Critical_buses=buslist(Critical_numbers);
	
	%Get rid of Xfrmr
	tic
	fprintf('\nRemoving unnecessary transformers:')
	
	%get all buses
	trf_buses=circuit.transformer(:).Buses;
	Store=[];
	xfrmrm_delete=[];
	for ii=1:length(trf_buses)
		bus_name=regexp(trf_buses{ii},'\.','split');
		
		%get indices of bus
		for jj=1:length(bus_name)
			buses=find(ismember(lower(buslist),lower(bus_name{jj}(1))));
			trf_bus_ind(ii,jj)=buses(1);
		end
		
		DownStream=[];
		for jj=2:size(find(trf_bus_ind(ii,:)>0),2)
			DownStream=[DownStream;cell2mat(generation(trf_bus_ind(ii,jj),3));trf_bus_ind(ii,jj)];
		end
		DownStream=unique(DownStream);
		
		CB=find(ismember(Critical_numbers,DownStream));
		if isempty(CB)
			
			%Bus on primary side
			Bus_num_prim=trf_bus_ind(ii,1);
			%nodes of primary bus
			Node_nums_prim=find(ismember(lower(YbusOrderVect),buslist(Bus_num_prim)));
			
			%Get nodes connected
			if circuit.transformer(ii).Phases==3
				Conn_nodes_prim=Node_nums_prim;
			else
				[~, keep]=strtok(trf_buses{ii},'\.'); conn_prim=keep(1);
				conn_prim=regexprep(conn_prim,'\.','','once'); conn_prim=regexp(conn_prim,'\.','split');
				Conn_nodes_prim=Node_nums_prim(find(ismember(YbusPhaseVect(Node_nums_prim),str2num(cell2mat([conn_prim{:}]')))));
			end
			
			for jj=2:size(trf_bus_ind,2)
				%nodes of secondary side
				Bus_num_sec=trf_bus_ind(ii,jj);
				%nodes of secondary bus
				Node_nums_sec=find(ismember(lower(YbusOrderVect),buslist(Bus_num_sec)));
				
				%Get nodes connected
				if circuit.transformer(ii).Phases==3
					Conn_nodes_sec(jj,:)=Node_nums_sec;
				else
					[~, keep]=strtok(trf_buses{ii},'\.');conn_sec=keep(jj); conn_sec=regexprep(conn_sec,'.0','','once');
					conn_sec=regexprep(conn_sec,'\.','','once'); conn_sec=regexp(conn_sec,'\.','split');
					Conn_nodes_sec(jj,:)=Node_nums_sec(find(ismember(YbusPhaseVect(Node_nums_sec),str2num(cell2mat([conn_sec{:}]')))));
				end
			end
			if any(Conn_nodes_sec==0)
				Conn_nodes_sec(find(Conn_nodes_sec==0))=[];
			end
			Remove_Ind_node=find(ismember(lower(YbusOrderVect),lower(buslist(DownStream))));
			
			Crit_Ind_node=Conn_nodes_prim;
			if ~isequal(unique(YbusPhaseVect(Remove_Ind_node)),unique(YbusPhaseVect(Conn_nodes_prim)))
				if pvFlag
					PV(Crit_Ind_node)=PV(Crit_Ind_node)+sum(PV(Remove_Ind_node));
					Weights(Crit_Ind_node,:)=Weights(Crit_Ind_node,:)+sum(Weights(Remove_Ind_node,:));
				end
				
				if loadFlag
					LOAD(Crit_Ind_node)=LOAD(Crit_Ind_node)+sum(LOAD(Remove_Ind_node));
				end
				
				I_cap(Crit_Ind_node)=I_cap(Crit_Ind_node)+sum(I_cap(Remove_Ind_node));
				
				if battFlag
					Batt(Crit_Ind_node,:)=Batt(Crit_Ind_node,:)+sum(Batt(Remove_Ind_node,:));
				end
				
				if MeasFlag
					InjP(Crit_Ind_node,1)=InjP(Crit_Ind_node,1)+InjP(Remove_Ind_node,1);
					InjP(Crit_Ind_node,2)=sqrt(InjP(Crit_Ind_node,2)^2+sum(InjP(Remove_Ind_node,2)^2));
					
					InjQ(Crit_Ind_node,1)=InjQ(Crit_Ind_node,1)+InjQ(Remove_Ind_node,1);
					InjQ(Crit_Ind_node,2)=sqrt(InjQ(Crit_Ind_node,2)^2+sum(InjQ(Remove_Ind_node,2)^2));
					
				end
				xfrmrm_delete=[xfrmrm_delete;(ii)];
				Store=[Store; DownStream];
			else
				% Match reduct bus with phase
				for jjj=1:length(Remove_Ind_node)
					
					MatchInd=find(ismember(YbusPhaseVect(Crit_Ind_node),YbusPhaseVect(Remove_Ind_node(jjj))))
					if pvFlag
						PV(Crit_Ind_node(MatchInd))=PV(Crit_Ind_node(MatchInd))+sum(PV(Remove_Ind_node));
						Weights(Crit_Ind_node(MatchInd),:)=Weights(Crit_Ind_node(MatchInd),:)+Weights(Remove_Ind_node(jjj),:);
					end
					
					if loadFlag
						LOAD(Crit_Ind_node(MatchInd))=LOAD(Crit_Ind_node(MatchInd))+LOAD(Remove_Ind_node(jjj));
					end
					
					I_cap(Crit_Ind_node(MatchInd))=I_cap(Crit_Ind_node(MatchInd))+I_cap(Remove_Ind_node(jjj));
					
					if battFlag
						Batt(Crit_Ind_node(MatchInd),:)=Batt(Crit_Ind_node(MatchInd),:)+Batt(Remove_Ind_node(jjj),:);
					end
					
					if MeasFlag
						InjP(Crit_Ind_node(MatchInd),1)=InjP(Crit_Ind_node(MatchInd),1)+InjP(Remove_Ind_node(jjj),1);
						InjP(Crit_Ind_node(MatchInd),2)=sqrt(InjP(Crit_Ind_node(MatchInd),2)^2+sum(InjP(Remove_Ind_node(jjj),2)^2));
						
						InjQ(Crit_Ind_node(MatchInd),1)=InjQ(Crit_Ind_node(MatchInd),1)+InjQ(Remove_Ind_node(jjj),1);
						InjQ(Crit_Ind_node(MatchInd),2)=sqrt(InjQ(Crit_Ind_node(MatchInd),2)^2+sum(InjQ(Remove_Ind_node(jjj),2)^2));
					end
					
				end
				xfrmrm_delete=[xfrmrm_delete;(ii)];
				Store=[Store; DownStream];
			end
		end
	end
	circuit.transformer(xfrmrm_delete)=[];
	trf_bus_ind(xfrmrm_delete,:)=[];
	
	%update everything
	[InjP,InjQ,MeasV,LineP,LineQ,Ybus,C_MAT, Lines,buslist,I_cap,Batt,PV,Weights,LOAD,YbusOrderVect,YbusPhaseVect,Node_number,Critical_buses,Store,Critical_numbers]=Update_Matrices(Critical_buses,Store,YbusOrderVect,buslist,Ybus,C_MAT,Lines,1,I_cap,battFlag,Batt,pvFlag,PV,Weights,loadFlag,MeasFlag,LOAD,YbusPhaseVect,InjP,InjQ,MeasV,LineP,LineQ);
	
	t_=toc;
	fprintf('time elapsed %f',t_)
	
	
	tic
	fprintf('\nAdding Remaining Transformer Nodes to CB: ')
	buses=[];
	trf=circuit.transformer;
	for ii=1:length(trf)
		buses{ii}=trf(ii).Buses;
		Bus1name=regexp(buses{ii}(1),'\.','split');
		bus1{ii}=char(Bus1name{1}(1));
		Bus2name=regexp(buses{ii}(2),'\.','split');
		bus2{ii}=char(Bus2name{1}(1));
	end
	bus1=regexp(bus1,'\.','split');
	bus1=[bus1{:}]'; bus1(find(cellfun('length',bus1)==1))=[];
	bus2=regexp(bus2,'\.','split');
	bus2=[bus2{:}]'; bus2(find(cellfun('length',bus2)==1))=[];
	xfrmrBuses=[bus1 bus2];
	[bus1_tmp, ia1, ic1]=unique(bus1);
	[bus2_tmp, ia2, ic2]=unique(bus2);
	keep1=find(ismember(lower(buslist),lower(bus1)));
	keep1=keep1(ic1);
	keep2=find(ismember(lower(buslist),lower(bus2)));
	keep2=keep2(ic2);
	New_Node_Num=reshape([keep1, keep2],[],1);
	New_Node_Num=unique(New_Node_Num);
	Critical_buses=[Critical_buses; bus1; bus2];
	Critical_numbers=[Critical_numbers; New_Node_Num];
	
	
	clear bus1 bus2 buses
	
	for ii=1:length(trf)
		if isempty(regexp(trf(ii).Name,'T'))
			continue
		else
			Critical_numbers=[Critical_numbers;cell2mat(generation(trf_bus_ind(ii,1),3))];
		end
	end
	Critical_numbers=unique(Critical_numbers); %get rid of repeat connections
	Critical_buses=buslist(Critical_numbers);
	t_=toc;
	fprintf('time elapsed %f',t_)
	
	% Get new topo
	tic
	fprintf('\nGetting New topogrophy: ')
	
	topo=zeros(max(Node_number),4);
	generation{1,1}=[];clear generation;  generation{1,1}=1; generation{1,4}=[];generation{1,5}=0;
	parent=1;
	topo(parent,1)=parent;
	[topo,generation]=topology_detect_large(topo,generation,Ybus,parent,Node_number);
	topo_view=topo;
	topo_view(find(topo_view(:,1)==0)',:)=[];
	c_new=0;
	t_=toc;
	fprintf('time elapsed %f',t_)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First find the nodes the nodes that connect two or more critical nodes,
%this will now also become a critical node. This is done by looking at the
%grandparents of each critical node and seeing if there are common nodes
%between two critical nodes. The closest common node is the one that is
%kept
tic
fprintf('\nGetting topogrophical critical Nodes (connection): ')
nn=length(Critical_numbers);
New_Critical_numbers=[];
for k=1:nn
	for j=k+1:nn
		%Here we find which parentso of CN(k) are also parents of CN(j),
		%return logic, find those indicices and return the parents of CN(k)
		VectorOfCommonParents=(generation{Critical_numbers(k),4}(find(ismember(cell2mat(generation(Critical_numbers(k),4)),cell2mat(generation(Critical_numbers(j),4))))));
		%find the node with the greatest distance (i.e. closest to each
		%node)
		[~,NewCriticalNumber]=max(cell2mat(generation(VectorOfCommonParents,5))); %Find the closest common point based on distance to substation (i.e. generation(_,5)
		New_Critical_numbers=[New_Critical_numbers, VectorOfCommonParents(NewCriticalNumber)];
	end
end
Critical_numbers=[Critical_numbers; New_Critical_numbers'];

%Add sourcebus to CN
Critical_numbers=[Critical_numbers; find(ismember(lower(buslist),'sourcebus'))];

%Add all substation equipment to list of critical buses
for ii=1:length(buslist)
	if ~(isempty(cell2mat(regexp(buslist(ii),'sub'))))
		Critical_numbers=[Critical_numbers; ii];
	end
end

Critical_numbers=unique(Critical_numbers); %get rid of repeat connections
Critical_buses=buslist(Critical_numbers);


% % % % treeplot(topo(:,2)')
% % % % [x,y] = treelayout(topo(:,2)');
% % % % for ii=1:length(x)
% % % %     text(x(ii),y(ii),num2str(ii))
% % % % end
t_=toc;
fprintf('time elapsed %f',t_)

if length(Critical_buses)~=length(Critical_numbers)
	stop=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Anything extra
% %Get line names to store lengths
% lineBus1=regexp(circuit.line{ii}.bus1,'\.','split','once');
% lineBus2=regexp(circuit.line{ii}.bus2,'\.','split','once');
% for ii=1:length(circuit.line)
% 	LineBus1{ii}=lineBus1{ii,1}(1);
% 	LineBus2{ii}=lineBus2{ii,1}(1);
% end

%Store Variables that are usefull later
AllCN=Critical_buses;
OriginalYbusOrderVect=YbusOrderVect;
OriginalYbusPhaseVect=YbusPhaseVect;
Originalbuslist=buslist;

fprintf('\n\nTotal critical nodes= %d\n',length(Critical_numbers))
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Start the reduction process w/ topology, critical nodes, and ybus
%Main part of code
%The algorithm is as follows:

%1) start at single critical node (CN), Check its children
%  a) if none are CN
%    i) delete all childrens rows of YBUS, combine all loads and PV to CN
%  b) if CN exist in children, repeat 1)
%2) start at single CN, check its parent
%  a) if its parent has >1 child
%    i) delete YBUS of children and collapse PV/load to parent
%  b) if parent has 1 child, go to next parent
%  c) repeat until parent is CN, go to 2)
%3) start at single CN,
%  a) check it's parent
%    i) if not CN
%      o) collapse parent with CN and its parent
%         o) Ybus adds, PV/Load weighted
%      o) move to next parent
%    ii) if CN
%      o) move to next CN
%Done?
Store=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 1) kill end nodes
%1) start at single critical node (CN), Check its children
%  a) if none are CN
%    i) delete all childrens rows of YBUS, combine all loads and PV to CN
%  b) if CN exist in children, repeat 1)
tic
fprintf('\nCollapsing End Nodes: ')

%Loop through al CN
for j=1:length(Critical_numbers)
	
	%Get the indices of the buses that are in need of deletion/aggregation
	EndBuses=cell2mat(generation(Critical_numbers(j),3));
	
	%Check to make sure that there are no cricitical nodes in it children
	if ~any(ismember(Critical_numbers,EndBuses));
		
		%It may be an end bus already, so skip
		if isempty(EndBuses)
			continue
		end
		
		%Get the node indices corresponding to the critical bus
		Ind=find(ismember(lower(YbusOrderVect),buslist(Critical_numbers(j))));
		
		%Alocate space to update load and node, actually not necessary due
		%to weighting
		if pvFlag
			PvMat=zeros(3,1);
		else
			pvMat=[];
		end
		if loadFlag
			LoadMat=zeros(3,1);
		else
			loadFlag=[];
		end
		I_mat=zeros(3,1);
		if battFlag
			batt_mat=zeros(3,4);
		else
			batt_mat=[];
		end
		if MeasFlag
			InjP_Mat=zeros(3,2);
			InjQ_Mat=zeros(3,2);
		else
			InjP_Mat=[];
			InjQ_Mat=[];
		end
		%Loop through buses that need to leave town
		for jj=1:length(EndBuses)
			
			% Find the Node indices corresponding to the end bus
			EndNodeInd=find(ismember(lower(YbusOrderVect),buslist(EndBuses(jj))));												% Find Indices of the individual end node
			
			% Find the indices of the phases that correspond to the end node and the
			% critical node
			EndPhasesInd=find(ismember(YbusPhaseVect(EndNodeInd),YbusPhaseVect(Ind)));												% Find the phases of that line, to keep only proper phases
			
			if pvFlag
				%Find the end nodes which have PV on them
				PvPhases=find(PV(EndNodeInd));
				
				%Match the CN Phases to update with the phases on end node that
				%have PV
				MatchPVandCN=find(ismember(YbusPhaseVect(Ind),YbusPhaseVect(EndNodeInd(PvPhases))));
				
				%Update weights
				Weights(Ind(MatchPVandCN),:)=Weights(Ind(MatchPVandCN),:)+Weights(EndNodeInd(PvPhases),:);																	% Update weights
				
				PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+PV(EndNodeInd(EndPhasesInd));			% Update Pv
				
			end
			if loadFlag
				LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+LOAD(EndNodeInd(EndPhasesInd));	% update load
			end
			I_mat(YbusPhaseVect(EndNodeInd))=I_mat(YbusPhaseVect(EndNodeInd))+I_cap(EndNodeInd);
			if battFlag
				batt_mat(YbusPhaseVect(EndNodeInd),:)=batt_mat(YbusPhaseVect(EndNodeInd),:)+batt_mat(EndNodeInd,:);
			end
			
			if MeasFlag
				InjP_Mat(YbusPhaseVect(EndNodeInd),1)=InjP_Mat(YbusPhaseVect(EndNodeInd),1)+InjP(EndNodeInd,1);
				InjP_Mat(YbusPhaseVect(EndNodeInd),2)=sqrt(InjP_Mat(YbusPhaseVect(EndNodeInd),2).^2+InjP(EndNodeInd,2).^2);
				InjQ_Mat(YbusPhaseVect(EndNodeInd),1)=InjQ_Mat(YbusPhaseVect(EndNodeInd),1)+InjQ(EndNodeInd,1);
				InjQ_Mat(YbusPhaseVect(EndNodeInd),2)=sqrt(InjQ_Mat(YbusPhaseVect(EndNodeInd),2).^2+InjQ(EndNodeInd,2).^2);
				
			end
			
		end
		
		if loadFlag
			LOAD(Ind)=LOAD(Ind)+LoadMat(YbusPhaseVect(Ind));
		end
		if pvFlag
			PV(Ind)=PV(Ind)+PvMat(YbusPhaseVect(Ind));
		end
		I_cap(Ind)=I_cap(Ind)+I_mat(YbusPhaseVect(Ind));
		if battFlag
			Batt(Ind,:)=Batt(Ind,:)+batt_mat(YbusPhaseVect(Ind),:);
		end
		if MeasFlag
			InjP(Ind,1)=InjP(Ind,1)+InjP_Mat(YbusPhaseVect(Ind),1);
			InjP(Ind,2)=sqrt(InjP(Ind,2).^2+InjP_Mat(YbusPhaseVect(Ind),2).^2);
			InjQ(Ind,1)=InjQ(Ind,1)+InjQ_Mat(YbusPhaseVect(Ind),1);
			InjQ(Ind,2)=sqrt(InjQ(Ind,2).^2+InjQ_Mat(YbusPhaseVect(Ind),2).^2);
		end
		
		%store children to delete at end. If you delte in real tiem it will
		%mess with the numbering.
		Store=[Store; EndBuses];
		
	end
end
t_=toc;
fprintf('time elapsed %f',t_)

%update everything
[InjP,InjQ,MeasV,LineP,LineQ,Ybus,C_MAT, Lines,buslist,I_cap,Batt,PV,Weights,LOAD,YbusOrderVect,YbusPhaseVect,Node_number,Critical_buses,Store,Critical_numbers]=Update_Matrices(Critical_buses,Store,YbusOrderVect,buslist,Ybus,C_MAT,Lines,1,I_cap,battFlag,Batt,pvFlag,PV,Weights,loadFlag,MeasFlag,LOAD,YbusPhaseVect,InjP,InjQ,MeasV,LineP,LineQ);

% Get new topo
tic
fprintf('\nGetting New topogrophy: ')

topo=zeros(max(Node_number),4);
generation{1,1}=[];clear generation;  generation{1,1}=1; generation{1,4}=[];generation{1,5}=0;
parent=1;
topo(parent,1)=parent;
[topo,generation]=topology_detect_large(topo,generation,Ybus,parent,Node_number);
topo_view=topo;
topo_view(find(topo_view(:,1)==0)',:)=[];
c_new=0;
t_=toc;
% % % treeplot(topo(:,2)')
% % % [x,y] = treelayout(topo(:,2)');
% % % for ii=1:length(x)
% % %     text(x(ii),y(ii),num2str(ii))
% % % end

fprintf('time elapsed %f',t_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 2: kill branches
%2) start at single CN, check its parent
%  a) if its parent has >1 child
%	 i) Remove CN child from list
%    ii) delete YBUS of children and collapse PV/load to parent
%  b) if parent has 1 child, go to next parent
%  c) repeat until parent is CN, go to 2)
tic
fprintf('\nCollapsing side branches: ')
for j=1:length(Critical_numbers)
	
	%get Parent of CN node
	AllGen=generation{Critical_numbers(j),3};
	
	%remove the children that are CN
	CNchildren=Critical_numbers(find(ismember(Critical_numbers,AllGen)));
	GenToRem=CNchildren;
	
	%remove branches that are attached to critical buses that don't have
	%any CN's in them
	for kk=1:length(CNchildren)
		Index=generation{CNchildren(kk),5}-generation{Critical_numbers(j),5};
		if Index>1
			CorrChild=generation{CNchildren(kk),4}(Index-1);
			GenToRem=[GenToRem;CorrChild;generation{CorrChild,3}];
		end
	end
	GenToRem=unique(GenToRem);
	GenToRemIndex=find(ismember(AllGen,GenToRem));
	AllGen(GenToRemIndex)=[];
	
	%find the nodes corresponding to the critical bus
	Ind=find(ismember(lower(YbusOrderVect),buslist(Critical_numbers(j))));
	if pvFlag
		PvMat=zeros(3,1);
	else
		PvMat=[];
	end
	if loadFlag
		LoadMat=zeros(3,1);
	else
		LoadMat=[];
	end
	I_mat=zeros(3,1);
	if battFlag
		batt_mat=zeros(3,4);
	else
		batt_mat=[];
	end
	if MeasFlag
		InjP_Mat=zeros(3,2);
		InjQ_Mat=zeros(3,2);
	else
		InjP_Mat=[];
		InjQ_Mat=[];
	end
	%Keep only corect phases
	for jj=1:length(AllGen)
		
		%Find the nodes that correspond to the bus being removed
		EndNodeInd=find(ismember(lower(YbusOrderVect),lower(buslist(AllGen(jj)))));												% Find Indices of the individual end node
		
		%find the indices end node phases that are also CN phases
		EndPhasesInd=find(ismember(YbusPhaseVect(EndNodeInd),YbusPhaseVect(Ind)));
		
		if pvFlag
			%Find the end nodes which have PV on them
			PvPhases=find(PV(EndNodeInd));
			
			%Match the CN Phases to update with the phases on end node that
			%have PV
			MatchPVandCN=find(ismember(YbusPhaseVect(Ind),YbusPhaseVect(EndNodeInd(PvPhases))));
			
			%Update weights
			Weights(Ind(MatchPVandCN),:)=Weights(Ind(MatchPVandCN),:)+Weights(EndNodeInd(PvPhases),:);																	% Update weights
			
			PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+PV(EndNodeInd(EndPhasesInd));			% Update Pv
			
		end
		if loadFlag
			LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+LOAD(EndNodeInd(EndPhasesInd));	% update load
		end
		I_mat(YbusPhaseVect(EndNodeInd))=I_mat(YbusPhaseVect(EndNodeInd))+I_cap(EndNodeInd);
		if battFlag
			batt_mat(YbusPhaseVect(EndNodeInd),:)=batt_mat(YbusPhaseVect(EndNodeInd),:)+batt_mat(EndNodeInd,:);
		end
		
		if MeasFlag
			InjP_Mat(YbusPhaseVect(EndNodeInd),1)=InjP_Mat(YbusPhaseVect(EndNodeInd),1)+InjP(EndNodeInd,1);
			InjP_Mat(YbusPhaseVect(EndNodeInd),2)=sqrt(InjP_Mat(YbusPhaseVect(EndNodeInd),2).^2+InjP(EndNodeInd,2).^2);
			InjQ_Mat(YbusPhaseVect(EndNodeInd),1)=InjQ_Mat(YbusPhaseVect(EndNodeInd),1)+InjQ(EndNodeInd,1);
			InjQ_Mat(YbusPhaseVect(EndNodeInd),2)=sqrt(InjQ_Mat(YbusPhaseVect(EndNodeInd),2).^2+InjQ(EndNodeInd,2).^2);
			
		end
	end
	
	%Update PV and Load
	if loadFlag
		LOAD(Ind)=LOAD(Ind)+LoadMat(YbusPhaseVect(Ind));
	end
	if pvFlag
		PV(Ind)=PV(Ind)+PvMat(YbusPhaseVect(Ind));
	end
	I_cap(Ind)=I_cap(Ind)+I_mat(YbusPhaseVect(Ind));
	if battFlag
		Batt(Ind,:)=Batt(Ind,:)+batt_mat(YbusPhaseVect(Ind),:);
	end
	if MeasFlag
		InjP(Ind,1)=InjP(Ind,1)+InjP_Mat(YbusPhaseVect(Ind),1);
		InjP(Ind,2)=sqrt(InjP(Ind,2).^2+InjP_Mat(YbusPhaseVect(Ind),2).^2);
		InjQ(Ind,1)=InjQ(Ind,1)+InjQ_Mat(YbusPhaseVect(Ind),1);
		InjQ(Ind,2)=sqrt(InjQ(Ind,2).^2+InjQ_Mat(YbusPhaseVect(Ind),2).^2);
	end
	
	%Store rows to delete
	Store=[Store; AllGen];
	
	
	%Remove branches along the main line that do not have any critical
	%nodes in them
	%Vector of parent nodes
	ParentNodes=generation(Critical_numbers(j),4);
	ParentNodes=ParentNodes{:};
	ParentNodes=[Critical_numbers(j);ParentNodes];
	
	for ii=2:length(ParentNodes)
		
		%This section is simply to remove the CN from the list of the
		%parents children
		AllGen=generation{ParentNodes(ii),3};
		
		%remove the parent you just came from
		GenToRem=find(ismember(generation{ParentNodes(ii),3},[ParentNodes(ii-1);generation{ParentNodes(ii-1),3}]));
		AllGen(GenToRem)=[];
		
		%Check to see if parent is critical node, if it is move to next CN
		if ismember(ParentNodes(ii),Critical_numbers)
			break
			
			%check to see if it has more children than just the CN.
		elseif (length(AllGen)<1) %go to next parent
			continue
			
			%If it has more children, then since it is not a cricitical node,
			%it should not have any CN in its children (otherwise would be
			%topographic CN
		else
			
			%find the children that have PV or laod and add them...basically
			%need to map between buslist and Ybus order
			Ind=find(ismember(lower(YbusOrderVect),buslist(ParentNodes(ii))));
			if pvFlag
				PvMat=zeros(3,1);
			else
				PvMat=[];
			end
			if loadFlag
				LoadMat=zeros(3,1);
			else
				LoadMat=[];
			end
			I_mat=zeros(3,1);
			if battFlag
				batt_mat=zeros(3,4);
			else
				batt_mat=[];
			end
			if MeasFlag
				InjP_Mat=zeros(3,2);
				InjQ_Mat=zeros(3,2);
			else
				InjP_Mat=[];
				InjQ_Mat=[];
			end
			%Keep only corect phases
			for jj=1:length(AllGen)
				
				%Find the end node indices
				EndNodeInd=find(ismember(lower(YbusOrderVect),lower(buslist(AllGen(jj)))));												% Find Indices of the individual end node
				
				%Find the phases of end node which match the CN
				EndPhasesInd=find(ismember(YbusPhaseVect(EndNodeInd),YbusPhaseVect(Ind)));
				
				if pvFlag
					%Find the end nodes which have PV on them
					PvPhases=find(PV(EndNodeInd));
					
					%Match the CN Phases to update with the phases on end node that
					%have PV
					MatchPVandCN=find(ismember(YbusPhaseVect(Ind),YbusPhaseVect(EndNodeInd(PvPhases))));
					
					%Update weights
					Weights(Ind(MatchPVandCN),:)=Weights(Ind(MatchPVandCN),:)+Weights(EndNodeInd(PvPhases),:);																	% Update weights
					
					PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+PV(EndNodeInd(EndPhasesInd));			% Update Pv
					
				end
				if loadFlag
					LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+LOAD(EndNodeInd(EndPhasesInd));	% update load
				end
				I_mat(YbusPhaseVect(EndNodeInd))=I_mat(YbusPhaseVect(EndNodeInd))+I_cap(EndNodeInd);
				if battFlag
					batt_mat(YbusPhaseVect(EndNodeInd),:)=batt_mat(YbusPhaseVect(EndNodeInd),:)+batt_mat(EndNodeInd,:);
				end
				
				if MeasFlag
					InjP_Mat(YbusPhaseVect(EndNodeInd),1)=InjP_Mat(YbusPhaseVect(EndNodeInd),1)+InjP(EndNodeInd,1);
					InjP_Mat(YbusPhaseVect(EndNodeInd),2)=sqrt(InjP_Mat(YbusPhaseVect(EndNodeInd),2).^2+InjP(EndNodeInd,2).^2);
					InjQ_Mat(YbusPhaseVect(EndNodeInd),1)=InjQ_Mat(YbusPhaseVect(EndNodeInd),1)+InjQ(EndNodeInd,1);
					InjQ_Mat(YbusPhaseVect(EndNodeInd),2)=sqrt(InjQ_Mat(YbusPhaseVect(EndNodeInd),2).^2+InjQ(EndNodeInd,2).^2);
					
				end
			end
			%Update PV and Load
			if loadFlag
				LOAD(Ind)=LOAD(Ind)+LoadMat(YbusPhaseVect(Ind));
			end
			if pvFlag
				PV(Ind)=PV(Ind)+PvMat(YbusPhaseVect(Ind));
			end
			I_cap(Ind)=I_cap(Ind)+I_mat(YbusPhaseVect(Ind));
			if battFlag
				Batt(Ind,:)=Batt(Ind,:)+batt_mat(YbusPhaseVect(Ind),:);
			end
			if MeasFlag
				InjP(Ind,1)=InjP(Ind,1)+InjP_Mat(YbusPhaseVect(Ind),1);
				InjP(Ind,2)=InjP(Ind,2)+InjP_Mat(YbusPhaseVect(Ind),2);
				InjQ(Ind,1)=InjQ(Ind,1)+InjQ_Mat(YbusPhaseVect(Ind),1);
				InjQ(Ind,2)=InjQ(Ind,2)+InjQ_Mat(YbusPhaseVect(Ind),2);
			end
			%Store rows to delete
			Store=[Store; AllGen];
		end
	end
end

%update everything
[InjP,InjQ,MeasV,LineP,LineQ,Ybus,C_MAT, Lines,buslist,I_cap,Batt,PV,Weights,LOAD,YbusOrderVect,YbusPhaseVect,Node_number,Critical_buses,Store,Critical_numbers]=Update_Matrices(Critical_buses,Store,YbusOrderVect,buslist,Ybus,C_MAT,Lines,1,I_cap,battFlag,Batt,pvFlag,PV,Weights,loadFlag,MeasFlag,LOAD,YbusPhaseVect,InjP,InjQ,MeasV,LineP,LineQ);

Reduction(2)=length(Store);

t_=toc;
fprintf('time elapsed %f',t_)

% Get new topo
tic
fprintf('\nGetting New topogrophy: ')
topo=zeros(max(Node_number),4);
generation{1,1}=[];clear generation;  generation{1,1}=1; generation{1,4}=[];generation{1,5}=0;
parent=1;
topo(parent,1)=parent;
[topo,generation]=topology_detect_large(topo,generation,Ybus,parent,Node_number);
topo_view=topo;
topo_view(find(topo_view(:,1)==0)',:)=[];
c_new=0;
t_=toc;
% % % treeplot(topo(:,2)')
% % % [x,y] = treelayout(topo(:,2)');
% % % for ii=1:length(x)
% % %     text(x(ii),y(ii),num2str(ii))
% % % end
fprintf('time elapsed %f',t_)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Part 3: Collapse main line down
%3) start at single CN,
%  a) check it's parent
%    i) if not CN
%      o) collapse parent with CN and its parent
%         o) Ybus adds, PV/Load weighted
%      o) move to next parent
%    ii) if CN
%      o) move to next CN

tic
fprintf('\nCollapsing Main Line: ')
Ycomb=strcat(YbusOrderVect,'.', num2str(YbusPhaseVect));
for j=length(Critical_numbers):-1:2
	
	%get Parent of CN node
	ParentNodes=generation(Critical_numbers(j),4); 	ParentNodes=ParentNodes{:};
	ParentNodes=[Critical_numbers(j);ParentNodes];
	
	for ii=2:length(ParentNodes)
		
		%Check if Parent is CN
		if ismember(ParentNodes(ii),Critical_numbers)
			break
			
			%If not CN
		else
			%Here we are getting the Z between each node to collape between
			%3 nodes, node 1, 2, and 3 to make two equivalent nodes. Z1 is
			%the connection between the CN (node 1) and the closest remaining
			%parent (node 2). Z2 is the connection between the closest remaining
			%parent and its parent (node 3). Zeq is the equivalent
			%connection after removign node 2. We then update the Ybus to
			%reflect these connections. At the end of the loop we remove
			%all extra rows and columns from Ybus. We also collapse the
			%PV/load by a weighted average.
			
			%Get indices.
			%The ind of the nodes to be deleted (middle node)
			IndMid=find(strcmpi(buslist(ParentNodes(ii)),YbusOrderVect))';
			
			%The ind of the critical nodes
			IndCrit=find(strcmpi(buslist(ParentNodes(1)),YbusOrderVect))';
			
			%The ind of the parent node on the other side
			IndUpStream=find(strcmpi(buslist(ParentNodes(ii+1)),YbusOrderVect))';
			
			%Find out if any missed connections make a matrix of
			%connections work correspondign to that.
			
			%phases             a    b    c
			%         Crit      x    x    x
			%         Middle    x    x    x
			%         UpStream  x    x    x
			
			ConnMat=zeros(3);
			ConnMat(1,YbusPhaseVect(IndCrit))=IndCrit;  %set row crit nodes
			ConnMat(2,YbusPhaseVect(IndMid))=IndMid;		%set row midd nodes
			ConnMat(3,YbusPhaseVect(IndUpStream))=IndUpStream;		%set row next nodes
			
			%Here we are just looping through to figure out connections
			ColFlag=zeros(3,1)';
			
			%Looping through the columns (phases) of the matrix
			for iii=1:3
				
				%If all three phases are connected, do nothing
				if length(find(ConnMat(:,iii)>0))==3
					ColFlag(iii)=1;
					
					%if two phases are connected, determine wo which sided
				elseif length(find(ConnMat(:,iii)>0))==2
					
					%if critical and middle, push everything to critical
					if ConnMat(1,iii)>0 && ConnMat(2,iii)>0
						if pvFlag
							PV(IndCrit(iii))=PV(IndCrit(iii))+PV(IndMid(iii));
							Weights(IndCrit(iii),:)=Weights(IndCrit(iii),:)+Weights(IndMid(iii),:);
						end
						if loadFlag
							LOAD(IndCrit(iii),:)=LOAD(IndCrit(iii),:)+LOAD(IndMid(iii),:);
						end
						if battFlag
							Batt(IndCrit(iii),:)=Batt(IndCrit(iii),:)+Batt(IndMid(iii),:);
						end
						if MeasFlag
							InjP(IndCrit(iii),1)=InjP(IndCrit(iii),1)+InjP(IndMid(iii),1);
							InjQ(IndCrit(iii),1)=InjQ(IndCrit(iii),1)+InjQ(IndMid(iii),1);
							InjP(IndCrit(iii),2)=sqrt(InjP(IndCrit(iii),2).^2+InjP(IndMid(iii),2).^2);
							InjQ(IndCrit(iii),2)=sqrt(InjQ(IndCrit(iii),2).^2+InjQ(IndMid(iii),2).^2);
						%%%consider adding removed line as Pinj Qinj
						end
						
						I_cap(IndCrit(iii),:)=I_cap(IndCrit(iii),:)+I_cap(IndMid(iii),:);

						%if middle and next, push everything to next node
					elseif ConnMat(3,1)>0 && ConnMat(2,1)>0
						if pvFlag
							PV(IndUpStream(iii))=PV(IndUpStream(iii))+PV(IndMid(iii));
							Weights(IndUpStream(iii),:)=Weights(IndUpStream(iii),:)+Weights(IndMid(iii),:);
						end
						if loadFlag
							LOAD(IndUpStream(iii))=LOAD(IndUpStream(iii))+LOAD(IndMid(iii));
						end
						if battFlag
							Batt(IndUpStream(iii),:)=Batt(IndUpStream(iii),:)+Batt(IndMid(iii),:);
						end
						if MeasFlag
							InjP(IndUpStream(iii),1)=InjP(IndUpStream(iii),1)+InjP(IndMid(iii),1);
							InjQ(IndUpStream(iii),1)=InjQ(IndUpStream(iii),1)+InjQ(IndMid(iii),1);
							InjP(IndUpStream(iii),2)=sqrt(InjP(IndUpStream(iii),2).^2+InjP(IndMid(iii),2).^2);
							InjQ(IndUpStream(iii),2)=sqrt(InjQ(IndUpStream(iii),2).^2+InjQ(IndMid(iii),2).^2);
						%%%consider adding removed line as Pinj Qinj
						end
						
						I_cap(IndUpStream(iii))=I_cap(IndUpStream(iii))+I_cap(IndMid(iii));

					end
				end
			end
			
			IndCritmod=ConnMat(1,find(ColFlag));
			IndMidMod=ConnMat(2,find(ColFlag));
			IndUpStreamMod=ConnMat(3,find(ColFlag));
			
			Z1=inv(Ybus(IndMidMod,IndCritmod));
			Z2=inv(Ybus(IndMidMod,IndUpStreamMod));
			
			%get equivalent Z
			Zeq=Z1+Z2;
			
			%get equivalent C
			C_MAT(IndUpStreamMod,IndCritmod)=C_MAT(IndMidMod,IndCritmod)+C_MAT(IndUpStreamMod,IndMidMod);
			
			%rewrite ybus
			Ybus(IndCritmod,IndUpStreamMod)=inv(Zeq);
			Ybus(IndUpStreamMod,IndCritmod)=inv(Zeq);
			
			%calculate weightings
			M2=Z2*inv(Zeq);	M1=Z1*inv(Zeq);
			
			
			%Update PV and Load
			if loadFlag
				LOAD(IndCritmod)=LOAD(IndCritmod)+(M2*LOAD(IndMidMod));
				LOAD(IndUpStreamMod)=LOAD(IndUpStreamMod)+(M1*LOAD(IndMidMod));
			end
				I_cap(IndCritmod)=I_cap(IndCritmod)+(M2*I_cap(IndMidMod));
				I_cap(IndUpStreamMod)=I_cap(IndUpStreamMod)+(M1*I_cap(IndMidMod));
			if pvFlag
				PV(IndCritmod)=PV(IndCritmod)+(M2*PV(IndMidMod));
				PV(IndUpStreamMod)=PV(IndUpStreamMod)+(M1*PV(IndMidMod));
				Weights(IndCritmod,:)=Weights(IndCritmod,:)+(M2*Weights(IndMidMod,:));
				Weights(IndUpStreamMod,:)=Weights(IndUpStreamMod,:)+(M1*Weights(IndMidMod,:));
			end
			if battFlag
				Batt(IndCritmod,:)=Batt(IndCritmod,:)+(M2*Batt(IndMidMod,:));
				Batt(IndUpStreamMod,:)=Batt(IndUpStreamMod,:)+(M1*Batt(IndMidMod,:));
			end
			if MeasFlag
				InjP(IndCritmod,1)=InjP(IndCritmod,1)+(M2*InjP(IndMidMod,1));
				InjQ(IndCritmod,1)=InjQ(IndCritmod,1)+(M2*InjQ(IndMidMod,1));
				InjP(IndUpStreamMod,1)=InjP(IndUpStreamMod,1)+(M2*InjP(IndMidMod,1));
				InjQ(IndUpStreamMod,1)=InjQ(IndUpStreamMod,1)+(M2*InjQ(IndMidMod,1));
				
				InjP(IndCritmod,2)=sqrt(InjP(IndCritmod,2).^2+(M2*InjP(IndMidMod,2)).^2);
				InjQ(IndCritmod,2)=sqrt(InjQ(IndCritmod,2).^2+(M2*InjQ(IndMidMod,2)).^2);
				InjP(IndUpStreamMod,2)=sqrt(InjP(IndUpStreamMod,2).^2+(M2*InjP(IndMidMod,2)).^2);
				InjQ(IndUpStreamMod,2)=sqrt(InjQ(IndUpStreamMod,2).^2+(M2*InjQ(IndMidMod,2)).^2);
			
				for jj=1:length(IndCritmod)
				% aggregate Line measurements
				 %find line with measurements to and from Critical node
				P_Crit_To=find(ismember(LineP(:,1),Ycomb(IndCritmod(jj))));
				P_To_Crit=find(ismember(LineP(:,2),Ycomb(IndCritmod(jj))));

				P_To_Mid=find(ismember(LineP(:,2),Ycomb(IndMidMod(jj))));
				P_Mid_To=find(ismember(LineP(:,1),Ycomb(IndMidMod(jj))));
				
				P_To_Up=find(ismember(LineP(:,2),Ycomb(IndUpStreamMod(jj))));
				P_Up_To=find(ismember(LineP(:,1),Ycomb(IndUpStreamMod(jj))));

				Line_Crit_To_Mid=P_Crit_To(find(ismember(P_Crit_To,P_To_Mid)));
				Line_Mid_To_Crit=P_Mid_To(find(ismember(P_Mid_To,P_To_Crit)));

				Line_Up_To_Mid=P_Up_To(find(ismember(P_Up_To,P_To_Mid)));
				Line_Mid_To_Up=P_Mid_To(find(ismember(P_Mid_To,P_To_Up)));
				
				LineP(end+1,1)=LineP(Line_Crit_To_Mid,1); LineP(end,2)=LineP(Line_Mid_To_Up,2);
				LineP(end,3)={cell2mat(LineP(Line_Crit_To_Mid,3))+cell2mat(LineP(Line_Mid_To_Up,3))};
				LineP(end,4)={sqrt(cell2mat(LineP(Line_Crit_To_Mid,4)).^2+cell2mat(LineP(Line_Mid_To_Up,4)).^2)};
				
				LineP(end+1,1)=LineP(Line_Up_To_Mid,1); LineP(end,2)=LineP(Line_Mid_To_Crit,2);
				LineP(end,3)={cell2mat(LineP(Line_Up_To_Mid,3))+cell2mat(LineP(Line_Mid_To_Crit,3))};
				LineP(end,4)={sqrt(cell2mat(LineP(Line_Up_To_Mid,4)).^2+cell2mat(LineP(Line_Mid_To_Crit,4)).^2)};	
				
				Delete=sort([Line_Crit_To_Mid,Line_Mid_To_Crit,Line_Up_To_Mid,Line_Mid_To_Up],'descend');
				LineP(Delete,:)=[]; 
				
				Q_Crit_To=find(ismember(LineQ(:,1),Ycomb(IndCritmod(jj))));
				Q_To_Crit=find(ismember(LineQ(:,2),Ycomb(IndCritmod(jj))));

				Q_To_Mid=find(ismember(LineQ(:,2),Ycomb(IndMidMod(jj))));
				Q_Mid_To=find(ismember(LineQ(:,1),Ycomb(IndMidMod(jj))));
				
				Q_To_Up=find(ismember(LineQ(:,2),Ycomb(IndUpStreamMod(jj))));
				Q_Up_To=find(ismember(LineQ(:,1),Ycomb(IndUpStreamMod(jj))));

				Line_Crit_To_Mid=Q_Crit_To(find(ismember(Q_Crit_To,Q_To_Mid)));
				Line_Mid_To_Crit=Q_Mid_To(find(ismember(Q_Mid_To,Q_To_Crit)));

				Line_Up_To_Mid=Q_Up_To(find(ismember(Q_Up_To,Q_To_Mid)));
				Line_Mid_To_Up=Q_Mid_To(find(ismember(Q_Mid_To,Q_To_Up)));
				
				LineQ(end+1,1)=LineQ(Line_Crit_To_Mid,1); LineQ(end,2)=LineQ(Line_Mid_To_Up,2);
				LineQ(end,3)={cell2mat(LineQ(Line_Crit_To_Mid,3))+cell2mat(LineQ(Line_Mid_To_Up,3))};
				LineQ(end,4)={sqrt(cell2mat(LineQ(Line_Crit_To_Mid,4)).^2+cell2mat(LineQ(Line_Mid_To_Up,4)).^2)};
				
				LineQ(end+1,1)=LineQ(Line_Up_To_Mid,1); LineQ(end,2)=LineQ(Line_Mid_To_Crit,2);
				LineQ(end,3)={cell2mat(LineQ(Line_Up_To_Mid,3))+cell2mat(LineQ(Line_Mid_To_Crit,3))};
				LineQ(end,4)={sqrt(cell2mat(LineQ(Line_Up_To_Mid,4)).^2+cell2mat(LineQ(Line_Mid_To_Crit,4)).^2)};	

				Delete=sort([Line_Crit_To_Mid,Line_Mid_To_Crit,Line_Up_To_Mid,Line_Mid_To_Up],'descend');
				LineQ(Delete,:)=[]; 
				end
			end
			
		end
		
		%Store stuff to delete so that it is gone!
		Store=[Store;ParentNodes(ii)];
		
		%get line lengths
		%line from parent node to CN                                %line to parent node from its parent
		tmp1=cell2mat(Lines(find(ismember(lower(Lines(:,2)),buslist(ParentNodes(ii)))),1))+cell2mat(Lines(find(ismember(lower(Lines(:,3)),buslist(ParentNodes(ii)))),1));
		tmp2=Lines(find(ismember(lower(Lines(:,3)),buslist(ParentNodes(ii)))),2);
		tmp3=Lines(find(ismember(lower(Lines(:,2)),buslist(ParentNodes(ii)))),3);
		Lines{end+1,1}=tmp1;
		Lines{end,2}=char(tmp2);
		Lines{end,3}=char(tmp3);
		Lines(find(ismember(lower(Lines(:,2)),buslist(ParentNodes(ii)))),:)=[];
		Lines(find(ismember(lower(Lines(:,3)),buslist(ParentNodes(ii)))),:)=[];
		
		
	end
end

%update everything
[InjP,InjQ,MeasV,~,~,Ybus,C_MAT, Lines,buslist,I_cap,Batt,PV,Weights,LOAD,YbusOrderVect,YbusPhaseVect,Node_number,Critical_buses,Store,Critical_numbers]=Update_Matrices(Critical_buses,Store,YbusOrderVect,buslist,Ybus,C_MAT,Lines,[],I_cap,battFlag,Batt,pvFlag,PV,Weights,loadFlag,MeasFlag,LOAD,YbusPhaseVect,InjP,InjQ,MeasV,LineP,LineQ);
Reduction(3)=length(Store);

t_=toc;
fprintf('time elapsed %f',t_)
%% Get new topo
tic
fprintf('\nGetting New topogrophy: ')
topo=zeros(max(Node_number),4);
generation{1,1}=[];clear generation;  generation{1,1}=1; generation{1,4}=[];generation{1,5}=0;
parent=1;
topo(parent,1)=parent;
[topo,generation]=topology_detect_large(topo,generation,Ybus,parent,Node_number);
topo_view=topo;
topo_view(find(topo_view(:,1)==0)',:)=[];
c_new=0;
t_=toc;
% % % treeplot(topo(:,2)')
% % % [x,y] = treelayout(topo(:,2)');
% % % for ii=1:length(x)
% % %     text(x(ii),y(ii),num2str(ii))
% % % end
fprintf('time elapsed %f',t_)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Finally: rewrite circuit into struct
%Here we rewrite the circuit based on the remaining nodes
%We rewrite the PV by using the updated kW for the nodes only. The
%remaining values are the same as default. Same for load.

	trf_bus=strtok(reshape([circuit.transformer{:}.buses],2,[]),'\.'); trf_bus(1,:)=[];
	[~,trf_bus_ind]=ismember(lower(trf_bus),lower(buslist));
	[trf_bus_ind, IA]=unique(trf_bus_ind); trf_bus=trf_bus(IA);
	trf_kV=cell2mat(reshape([circuit.transformer{:}.Kv],2,[])); trf_kV(1,:)=[];
	trf_phase=[circuit.transformer{:}.Phases];
	trf_phase(find(trf_phase==3))=4; trf_phase(find(trf_phase==1))=3; trf_phase(find(trf_phase==4))=1;
	trf_kV=trf_kV.*sqrt(trf_phase);  trf_kV=trf_kV(IA);

tic
fprintf('\n\nFinally rewriting the circuit: ')
if pvFlag
	%PV
	%get PV kVa and pf
	NodeWithPV=find(PV>0);
	circuit.pvsystem=dsspvsystem;
	% BusListNumWithPV=find(ismember(buslist,unique(YbusOrderVect(BusWithPV))));
	for ii=1:length(NodeWithPV)
		circuit.pvsystem(ii)=dsspvsystem;
		% 	Ind=find(ismember(YbusOrderVect,buslist(BusListNumWithPV(ii))));
		% 	PhaseStr=[];
		% 	for j=1:length(Ind)
		% 		PhaseStr=[PhaseStr '.' num2str(YbusPhaseVect(Ind(j)))];
		% 	end
		circuit.pvsystem(ii).Name=['PV_on_' char(YbusOrderVect(NodeWithPV(ii))) '.' num2str(YbusPhaseVect(NodeWithPV(ii)))];%char(buslist(BusListNumWithPV(ii))) PhaseStr];
		circuit.pvsystem(ii).phases=1;
		circuit.pvsystem(ii).bus1=[char(YbusOrderVect(NodeWithPV(ii))) '.' num2str(YbusPhaseVect(NodeWithPV(ii)))];%[char(buslist(BusListNumWithPV(ii))) PhaseStr];
		circuit.pvsystem(ii).irradiance=1;
		circuit.pvsystem(ii).cutin=0;
		circuit.pvsystem(ii).cutout=0;
		circuit.pvsystem(ii).Pmpp=num2str(real(PV(NodeWithPV(ii))));
		circuit.pvsystem(ii).pf=num2str(real(PV(NodeWithPV(ii)))./sqrt(real(PV(NodeWithPV(ii))).^2+imag(PV(NodeWithPV(ii))).^2));
		circuit.pvsystem(ii).kVA=num2str(abs(PV(NodeWithPV(ii))));
				
		BusWithPV=find(ismember(buslist,lower(YbusOrderVect(NodeWithPV(ii)))));
		MatchInd=find(ismember(trf_bus_ind,[cell2mat(generation(BusWithPV,4)); BusWithPV]));
		circuit.pvsystem(ii).Kv=trf_kV(MatchInd(end))/sqrt(3);
	end
end
%% Load
%get PV kVa and pf
if loadFlag
	

	
	NodeWithLOAD=find(LOAD>0);
	circuit.load=dssload;
	% BusListNumWithLOAD=find(ismember(buslist,unique(YbusOrderVect(BusWithLOAD))));
	for ii=1:length(NodeWithLOAD)
		circuit.load(ii)=dssload;
		% 	Ind=find(ismember(YbusOrderVect,buslist(BusListNumWithLOAD(ii))));
		% 	PhaseStr=[];
		% 	for j=1:length(Ind)
		% 		PhaseStr=[PhaseStr '.' num2str(YbusPhaseVect(Ind(j)))];
		% 	end
		circuit.load(ii).Name=['LOAD_on_' char(YbusOrderVect(NodeWithLOAD(ii))) '.' num2str(YbusPhaseVect(NodeWithLOAD(ii)))];%['Load_on_' char(buslist(BusListNumWithLOAD(ii))) PhaseStr];
		circuit.load(ii).phases=1;
		circuit.load(ii).bus1=[char(YbusOrderVect(NodeWithLOAD(ii))) '.' num2str(YbusPhaseVect(NodeWithLOAD(ii)))];
		% 	circuit.load(ii).Pf=num2str(real(LOAD(BusWithLOAD(ii)))./sqrt(real(LOAD(BusWithLOAD(ii))).^2+imag(LOAD(BusWithLOAD(ii))).^2));
		circuit.load(ii).Kw=num2str(real(LOAD(NodeWithLOAD(ii))));
		circuit.load(ii).Kvar=num2str(imag(LOAD(NodeWithLOAD(ii))));
		circuit.load(ii).model=5;
		%Make sure kV is correct
		BusWithLoad=find(ismember(buslist,lower(YbusOrderVect(NodeWithLOAD(ii)))));
		MatchInd=find(ismember(trf_bus_ind,[cell2mat(generation(BusWithLoad,4)); BusWithLoad]));
		circuit.load(ii).Kv=trf_kV(MatchInd(end))/sqrt(3);
	end
	
	%Need to come up with elegant way of handling capacitive loads. Right
	%now, alll capacitive loads are being written to loads, but actrually
	%some of the capacitance is transferred through the cmatrix in step 3,
	%and thus is being double counted!!!!
% % 	%Load from capacitance
% % 	Loads_ind=find(abs(I_cap)>0);
% % 	
% % 	for iii=1:length(Loads_ind)
% % 		circuit.load(ii+iii)=dssload;
% % 		circuit.load(ii+iii).Name=['Capacitive_LOAD_on_' char(YbusOrderVect(Loads_ind(iii))) '.' num2str(YbusPhaseVect(Loads_ind(iii)))];
% % 		circuit.load(ii+iii).phases=1;
% % 		circuit.load(ii+iii).bus1=[char(YbusOrderVect(Loads_ind(iii))) '.' num2str(YbusPhaseVect(Loads_ind(iii)))];
% % 		circuit.load(ii+iii).Kw=0;
% % 		circuit.load(ii+iii).Kvar=num2str(I_cap(Loads_ind(iii)));
% % 		
% % 		%Make sure kV is correct
% % 		BusWithLoad=find(ismember(buslist,lower(YbusOrderVect(Loads_ind(iii)))));
% % 		MatchInd=find(ismember(trf_bus_ind,[cell2mat(generation(BusWithLoad,4)); BusWithLoad]));
% % 		circuit.load(ii+iii).Kv=trf_kV(MatchInd(end))/sqrt(3);
% % 	end
end
%% battery

if battFlag
	batts_ind=find(Batt(:,1)>0)
	
	for ii=1:length(batts_Ind)
		circuit.storage(ii)=dssstorage;
		circuit.storage(ii).Name=['storage'];
	end
end
%% delete useless classes now
names=fieldnames(circuit);
keepFields={'load','buslist','line','circuit','capcontrol','transformer','capacitor','basevoltages','regcontrol','pvsystem','reactor'};
names=names(find(~ismember(names,keepFields)));
for ii=1:length(names)
	circuit=rmfield(circuit,names{ii});
end

%Update circuit info
circuit.circuit.Name=[circuit.circuit.Name '_Reduced'];
tmp=find(ismember(lower(circuit.buslist.id),buslist));
circuit.buslist.id=circuit.buslist.id(tmp);
circuit.buslist.coord=circuit.buslist.coord(tmp,:);

count=0;
for ii=1:length(buslist)
	busTo=buslist(generation{ii,1});
	BusFrom=buslist(generation{ii,2});
	for jj=1:length(BusFrom)
		if length(BusFrom)>1
			stop=1;
		end
		count=count+1;
		Buses{count,1}=char(busTo);
		Buses{count,2}=char(BusFrom{jj});
	end
end
% Store=zeros(length(Buses),1);
% for ii=1:length(Buses)
% 	if all(ismember([Buses{ii,1} Buses{ii,2}],lower(xfrmrBuses)))
% 		Store(ii)=1;
% 	end
% end
% Buses(find(Store),:)=[];
% Store=[];
% for ii=1:length(Buses)
% 		Ind1=find(ismember(lower(LengthVect(:,1)),Buses{ii,1}));
% 		Ind2=find(ismember(lower(LengthVect(:,2)),Buses{ii,2}));
% 		LineLength(ii)=LengthVect(Ind1(find(ismember(Ind1,Ind2))),3);
% end
% LineLength=cell2mat(LineLength);
% Length(find([Length{:,2}]==0),:)=[];
% for ii=1:length(Buses)
% 	LineLength(find(ismember(Buses(:,2),[Length{ii,1}])))=cell2mat(Length(ii,2));
% end
trf_buses=circuit.transformer(:).Buses;
for ii=1:length(trf_buses)
	bus_name=regexp(trf_buses{ii},'\.','split');
	
	for jj=1:length(bus_name)
		buses=find(ismember(lower(buslist),lower(bus_name{jj}(1))));
		trf_bus_ind(ii,jj)=buses(1);
	end
	
	Bus_num=trf_bus_ind(ii,1);
	Node_nums=find(ismember(lower(YbusOrderVect),buslist(Bus_num)));
end

if isfield(circuit,'reactor')
	for ii=1:length(circuit.reactor)
		trf_bus_ind(end,1)=find(ismember(lower(buslist),lower(strtok(circuit.reactor{ii}.bus1,'.'))));
		trf_bus_ind(end,2)=find(ismember(lower(buslist),lower(strtok(circuit.reactor{ii}.bus2,'.'))));
	end
end


circuit.line=dssline;
count=0;
for ii=1:length(Buses)
	
	bus1Ind=find(ismember(buslist,Buses{ii,1}));
	bus2Ind=find(ismember(buslist,Buses{ii,2}));
	
	Match1=find(ismember(trf_bus_ind(:,1),bus1Ind));
	Match2=find(ismember(trf_bus_ind(:,2),bus2Ind));
	
	if ~isempty(find(ismember(Match1,Match2)))
		continue
	else
		count=count+1;
		circuit.line(count)=dssline;
		circuit.line(count).Name=[char(Buses{ii,1}) '_' char(Buses{ii,2})];
% 		circuit.line(count).Linecode=circuit.linecode(count).Name;
		circuit.line(count).Units='kft';
		circuit.line(count).R1=[];
		circuit.line(count).R0=[];
		circuit.line(count).X0=[];
		circuit.line(count).X1=[];
		
		Phases1=YbusPhaseVect(find(ismember(lower(YbusOrderVect),lower(Buses{ii,1}))));
		Phases2=YbusPhaseVect(find(ismember(lower(YbusOrderVect),lower(Buses{ii,2}))));
		Phases=Phases1(find(ismember(Phases1,Phases2)));
		PhaseStr=[];
		for j=1:length(Phases)
			PhaseStr=[PhaseStr '.' num2str(Phases(j))];
		end
		
		circuit.line(count).bus1=[char(Buses{ii,1}) PhaseStr];
		circuit.line(count).Phases=length(Phases);
		circuit.line(count).bus2=[char(Buses{ii,2}) PhaseStr];
		id1=find(ismember(lower(Lines(:,2)),char(Buses{ii,1})));
		id2=find(ismember(lower(Lines(:,3)),char(Buses{ii,2})));
		if length(find(ismember(id1,id2)))>1
			circuit.line(count).Length=Lines(id1(1),1);
		else
			circuit.line(count).Length=Lines(id1(find(ismember(id1,id2))),1);
		end
		
		%get r and x matrix
		Ind1= find(ismember(lower(YbusOrderVect),Buses{ii,1}));
		Ind2=find(ismember(lower(YbusOrderVect),Buses{ii,2}));
		
		%Make sure Phases are in correct order
		[~,I]=sort(YbusPhaseVect(Ind1));
		Ind1=Ind1(I);
		[~,I]=sort(YbusPhaseVect(Ind2));
		Ind2=Ind2(I);
		
		Mat=full(Ybus(Ind1,Ind2));
		missingrow=~any(Mat,2); missingcol=~any(Mat,1);
		Mat(missingrow,:)=[]; Mat(:,missingcol)=[];
		FullMat=-inv(Mat);
		
		id1=find(ismember(lower(Lines(:,2)),char(Buses{ii,1})));
		id2=find(ismember(lower(Lines(:,3)),char(Buses{ii,2})));
		
		if isempty(id1)
			FullMat=FullMat;
		elseif length(find(ismember(id1,id2)))>1
			FullMat=FullMat./cell2mat(Lines(id1(1),1));
			Line_Cap=C_MAT(Ind1,Ind2)./cell2mat(Lines(id1(1),1));
		else
			FullMat=FullMat./cell2mat(Lines(id1(find(ismember(id1,id2))),1));
			Line_Cap=C_MAT(Ind1,Ind2)./cell2mat(Lines(id1(find(ismember(id1,id2))),1));
		end
		
		
		if length(FullMat)==1
			circuit.line(count).Rmatrix=['(' num2str(real(FullMat(1,1))) ')'];
			circuit.line(count).Xmatrix=['(' num2str(imag(FullMat(1,1))) ')'];
			circuit.line(count).Cmatrix=['(' num2str(Line_Cap(1,1)) ')'];
		elseif length(FullMat)==2
			circuit.line(count).Rmatrix=['(' num2str(real(FullMat(1,1))) '|' num2str(real(FullMat(2,1:2)))  ')'];
			circuit.line(count).Xmatrix=['(' num2str(imag(FullMat(1,1))) '|' num2str(imag(FullMat(2,1:2)))  ')'];
			circuit.line(count).Cmatrix=['(' num2str(Line_Cap(1,1)) '|' num2str(Line_Cap(2,1:2)) ')'];
		else
			circuit.line(count).Rmatrix=['(' num2str(real(FullMat(1,1))) '|' num2str(real(FullMat(2,1:2))) '|' num2str(real(FullMat(3,1:3))) ')'];
			circuit.line(count).Xmatrix=['(' num2str(imag(FullMat(1,1))) '|' num2str(imag(FullMat(2,1:2))) '|' num2str(imag(FullMat(3,1:3))) ')'];
			circuit.line(count).Cmatrix=['(' num2str(Line_Cap(1,1)) '|' num2str(Line_Cap(2,1:2)) '|' num2str(Line_Cap(3,1:3)) ')'];
		end
	end

end

% % %update capcontrol to match updated line names
% % if isfield(circuit,'capcontrol')
% % 	for ii=1:length(capConBuses1)
% % 		Numms1=find(ismember(lower({Buses{:,1}}),lower(char(capConBuses1{ii}))));
% % 		Numms2=find(ismember(lower(Buses(:,2)),lower(char(capConBuses2{ii}))));
% % 		circuit.capcontrol(ii).Element=['line.' circuit.line(Numms1(find(ismember(Numms1,Numms2)))).Name];
% % 	end
% % end

%% 
Ycomb=strcat(YbusOrderVect,'.', num2str(YbusPhaseVect));
if MeasFlag
	n=length(InjP);
	m=length(InjQ);
	k=length(LineP);
	l=length(LineQ);
	j=length(MeasV);
	
	%Get real and reactive to correct parts
	Real_Q=real(InjQ(:,1));Imag_Q=imag(InjQ(:,1));
	Real_P=real(InjP(:,1));Imag_P=imag(InjP(:,1));
	
	InjQ(:,1)=Real_Q+Imag_P;
	InjP(:,1)=Real_P+Imag_Q;
	
	%Consider changing
	InjQ(:,2)=real(InjQ(:,2));
	InjP(:,2)=real(InjP(:,2));
	
	Measurements=num2cell(zeros(n+m+k+l+j,5));
	
	Measurements(1:n,1)={1};
	Measurements(1:n,2)=Ycomb; Measurements(1:n,3)={0};
	Measurements(1:n,4)=num2cell(InjP(:,1));
	Measurements(1:n,5)=num2cell(InjP(:,2));
	
	Measurements(n+1:n+m,1)={2};
	Measurements(n+1:n+m,2)=Ycomb; Measurements(n+1:n+m,3)={0};
	Measurements(n+1:n+m,4)=num2cell(InjQ(:,1));
	Measurements(n+1:n+m,5)=num2cell(InjQ(:,2));
	
	Measurements(n+1+m:n+m+k,1)={3};
	Measurements(n+1+m:n+m+k,2:5)=LineP;
	
	Measurements(n+1+m+k:n+m+k+l,1)={3};
	Measurements(n+1+m+k:n+m+k+l,2:5)=LineQ;
	
	Measurements(n+1+m+k+l:n+m+k+l+j,1)={5};
	Measurements(n+1+m+k+l:n+m+k+l+j,2)=Ycomb;
	Measurements(n+1+m+k+l:n+m+k+l+j,3)={0};
	Measurements(n+1+m+k+l:n+m+k+l+j,4)=num2cell(MeasV(:,1));
	Measurements(n+1+m+k+l:n+m+k+l+j,5)=num2cell(MeasV(:,2));
end

if MeasFlag
circuit.Measurements=Measurements;
end

%write weights
if pvFlag
	circuit.weights=Weights;
	circuit.PvbusMap=PVbusMap;
	circuit.PvSize=PvSize;
end
circuit.Ybus=full(Ybus);
circuit.Zbus=inv(full(Ybus));
circuit.pv=PV;
circuit.ld=LOAD;
circuit.CriticalNode=Ycomb;
% circuit.PhaseReduction=Reduction;

%% End
t_=toc;
fprintf('time elapsed %f \n',t_)

Ld_a_end=sum(LOAD(find(YbusPhaseVect==1)));
Ld_b_end=sum(LOAD(find(YbusPhaseVect==2)));
Ld_c_end=sum(LOAD(find(YbusPhaseVect==3)));

outputdss=WriteDSS(circuit,[],0,pwd);
stop_bufer=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [topo, generation]=topology_detect_large(topo,generation,Ybus,parent,Node_number)
% Function to get the feeder topology
% Created by Vahid R. Disfani

% topo:
% column 1: node #
% column 2: parent node #
% column 3: Number of children
% column 4: Number of downstream buses (grand children and so on)

% generation:
% column 1: Node #
% column 2: List of children
% column 3: List of downstream buses (grandchildren to the end points)
% column 4: List of grandparent nodes until the substation
% column 5: Distance from substation assuming that the distance of each line is 1 (number of grandparent nodes)

[Nbus,Nbuscol]=size(Ybus);
%get the nodes corresponding to the parent nodes
nodes=find(Node_number==parent);
%find nodes which are connected to parent nodes...basically need to do this
%to detect phases. Then keep only the uniqe nodes (i.e the ones that are
%connected)
adj_nodes=mod(find(Ybus(:,nodes)~=0)-.5,Nbus)+.5;
adj_bus=Node_number(adj_nodes);
[b1,m1,n1]=unique(adj_bus,'first');

adj_bus=adj_bus(sort(m1));
adj_bus(find(adj_bus==parent))=0; %delete parent from list
for i=1:max(Node_number)
	adj_bus(find(adj_bus==topo(i,2)))=0;
end
adj_bus(find(adj_bus==0))=[]; %remove parent
generation{parent,2}=[];	  %set children empty
generation{parent,3}=[];	  %set grandchildren empty
if length(adj_bus~=0)
	for k=1:length(adj_bus)
		child=adj_bus(k);
		if max(max(topo(:,1:2)==child))==0
			topo(child,1)=child;
			topo(child,2)=parent;
			generation{child,1}=child;
			generation{child,4}=[parent;generation{parent,4}];
			generation{child,5}= generation{parent,5}+1;
			[topo, generation]=topology_detect_large(topo,generation,Ybus,child,Node_number);
			topo(parent,3)=topo(parent,3)+1;
			topo(parent,4)=topo(parent,4)+topo(child,4)+1;
			generation{parent,2}=[generation{parent,2};child];
			generation{parent,3}=[generation{parent,3};child;generation{child,3}];
		end
	end
end
end

function [InjP,InjQ,MeasV,LineP,LineQ,Ybus,C_MAT, Lines,buslist,I_cap,Batt,PV,Weights,LOAD,YbusOrderVect,YbusPhaseVect,Node_number,Critical_buses,Store,Critical_numbers]=Update_Matrices(Critical_buses,Store,YbusOrderVect,buslist,Ybus,C_MAT,Lines,LinesUpdate,I_cap,battFlag,Batt,pvFlag,PV,Weights,loadFlag,MeasFlag,LOAD,YbusPhaseVect,InjP,InjQ,MeasV,LineP,LineQ)

%Update store to reflect Ybus order
Store=unique(Store);
S=find(ismember(lower(YbusOrderVect),buslist(Store)));

%Reduceeverything to collapse chitlins
Ybus(S,:)=[]; Ybus(:,S)=[];
C_MAT(S,:)=[]; C_MAT(:,S)=[];

if ~isempty(LinesUpdate)
	Lines(find(ismember(Lines(:,2),buslist(Store))),:)=[];
	Lines(find(ismember(Lines(:,3),buslist(Store))),:)=[];
end

if MeasFlag
%Reduce everything to collapse chitlins
InjP(S,:)=[]; InjQ(S,:)=[]; MeasV(S,:)=[];

LineP(find(ismember(strtok(LineP(:,1),'.'),buslist(Store))),:)=[];
LineP(find(ismember(strtok(LineP(:,2),'.'),buslist(Store))),:)=[];
LineQ(find(ismember(strtok(LineQ(:,1),'.'),buslist(Store))),:)=[];
LineQ(find(ismember(strtok(LineQ(:,2),'.'),buslist(Store))),:)=[];
end
I_cap(S)=[];
if battFlag
	Batt(S,:)=[];
end
if pvFlag
	PV(S)=[];
	Weights(S,:)=[];
end
if loadFlag
	LOAD(S,:)=[];
end
YbusOrderVect(S)=[];
YbusPhaseVect(S)=[];
Node_number=[];
buslist(Store)=[];
Store=[];

for ii=1:length(buslist)
	Ind=find(strcmpi(buslist(ii),YbusOrderVect))';
	Node_number(Ind)=ii;
end

%Need to update critical numbers to busnames
Critical_numbers=find(ismember(buslist,Critical_buses));
end

function [Output]=Convert_to_kft(units, Input)

if strcmp(units,'kft')
	Output=Input;
elseif strcmp(units,'mi')
	Output=Input*5.280;
elseif strcmp(units,'km')
	Output=Input*3.28084;
elseif strcmp(units,'m')
	Output=Input*.00328084;
elseif strcmp(units,'ft')
	Output=Input/1000;
elseif strcmp(units,'in')
	Output=Input/12/1000;
elseif strcmp(units,'cm')
	Output=Input/2.54/12/1000;
elseif strcmp(units,'none')
	Output=Input;
end
end

function [cir cmds] = dssparse(filename)
% Parse OpenDSS file to OpenDSS struct in Matlab
% Outputs:
%			cir : openDSS circuit struct with all components
%			cmds : list of commands/settings to run simulation

%process inputs
if ischar( filename )
	id = find((filename=='/')|(filename=='\'),1,'last');
	fdir = filename(1:id);
	if(~strcmp(filename(1:2),'\\') && filename(2)~=':')
		fdir = [pwd '/' fdir];
		filename = [pwd '/' filename];
	end
else
	error('Invalid input. Must specify filename.')
end

% initialize
cir = struct();
cmds = '';
warningnames = {'dsscapacitor:grounding','cleanPhase:phaselargerthan3'};
for i=1:length(warningnames)
	oldwarnings(i) = warning('off',warningnames{i});
end


% load file
fid = fopen(filename);

ignored = 0;
knownobjs = {};
unknownobjs = {};

while 1
    try
        l = fgetl(fid);
    catch e
        error([filename ' doesn''t exist!']);
    end
    if ~ischar(l),   break,   end
    
	% remove all spaces/tabs at the beginning of the line
	l = strtrim(l);
	
	% skip comments and empty lines; remove any trailing comments
	l = regexp(l,'!|//','split');
	% clean up spaces
	l = regexprep(l,'\s*=\s*','=');
	l = regexprep(l,'\s+',' ');
	
	l = l{1};
	if strcmp(l,'')
		continue;
	end
	% remove 'object=' string if exist
    if ~isempty(strfind(lower(l),'object='))
        i = strfind(lower(l),'object=');
        l(i:i+6) = [];
    end
    
	% Handle the meat content!!!
	cmd = regexp(l,'\S+','match','once');
            
	switch lower(cmd)
		% ignore 'clear' command
		case 'clear'
			continue;
		% handle 'new' command
		case 'new'
            % search for class name and object name
			n  = regexpi(l,'(\S+)\.(\S+)','once','tokens');
			cn = n{1};
            
			if ~ismember(cn,knownobjs) && ~ismember(cn,unknownobjs)
				try
					feval(['dss' lower(cn)]);
					knownobjs = [knownobjs cn];
				catch err
					unknownobjs = [unknownobjs cn];
				end
			end
			
			if ismember(cn,knownobjs)
				ignored = 0;
				obj = createObj(l);
				cn = class(obj);
				cn = cn(4:end);
				if ~isfield(cir,cn)
					cir.(cn) = obj;
				else
					cir.(cn)(end+1) = obj;
				end
			else
				ignored = 1;
			end
		% handle lines start with '~' (continuing of the previous "new" command)
		case '~'
			if ~ignored 
				obj = addtoObj(obj,l);
				cn = class(obj);
				cn = cn(4:end);
				cir.(cn)(end) = obj;
			end
		% handle 'set' command
		case 'set'
			% special handle for base voltages
			if strfind(lower(l),'voltagebases')
				basev = regexp(lower(l),'voltagebases=([\[\(\{"''][^=]+[\]\)\}''"]|[^"''\[\(\{]\S*)','tokens');
				basev = regexp(basev{1}{1},'[\d.]+','match');
				cir.basevoltages = cellfun(@str2num,basev);
			else
				cmds = sprintf('%s\n%s',cmds,l);
			end
		% handle 'redirect' commands
		case {'redirect','compile'}
			fn = regexp(l,'\s+','split');
			fn = strtrim(fn);
			if length(fn) < 2
				error('Check redirect/compile command'); 
			else
				fn = fn{2};
			end
			% check if file name is wrapped in quotes
			m = regexp(l,'[\(\[\{"\''](.*)[\)\]\}"\'']','tokens');
			
			if ~isempty(m)
				fn = m{1}{1};
			end
			
			% get subcir and sub command list
			[cir2 cmdlist] = dssparse([fdir '/' fn]);
			% merge sub circuit to original circuit
			cmds = sprintf('%s\n%s',cmds,cmdlist);
			fnames = fieldnames(cir2);
			for i = 1:length(fnames)
				fn_ = fnames{i};
				if ismember(fnames(i),fieldnames(cir))
					cir.(fn_) = [cir.(fn_) cir2.(fn_)];
				else
					cir.(fn_) = cir2.(fn_);
				end
			end
		% handle other commands
		case 'buscoords' 
			bfn = regexp(l,' ','split');
			bfn = bfn{2};
			fid2 = fopen([fdir '/' bfn]);
			try 
				dat = textscan(fid2,'%[^,]%*[,]%f%*[,]%f%*[\r\n]');
			catch err
				try 
					fseek(fid2,0,'bof');
					dat = textscan(fid2,'%s%f%f%*[\r\n ]');
				catch err
					error('dssparse:buscoordsfileinput','not recoganized input file format. Please use comma or space as delimiter');
				end
			end
			cir.buslist.id = dat{:,1};
			cir.buslist.coord = [dat{:,2}, dat{:,3}];
			fclose(fid2);
		otherwise
			% add to cmds
			cmds = sprintf('%s\n%s\n',cmds,l);
	end
end

for i=1:length(oldwarnings)
	warning(oldwarnings(i).state,oldwarnings(i).identifier);
end

try
	cir.switch = cir.swtcontrol;
	cir = rmfield(cir,'swtcontrol');
catch
end

% close file
fclose(fid);
if ~isempty(unknownobjs)
	disp('Unimplemented object(s):');
	disp(unknownobjs);
end
end

function obj = createObj(l)

% search for class name and object name
[n,prop]  = regexpi(l,'(\S+)\.(\S+)','once','tokens','split');
cn = n{1};
on = n{2};

% create object
obj = feval(['dss' lower(cn)]);
obj.Name = on;

% add properties
if ~isempty(prop{2})
	obj = addtoObj(obj,prop{2});
end

end

function obj = addtoObj(obj,l)

props = regexp(l,'(\S+)=([\[\(\{"''][^=]+[\]\)\}''"]|[^"''\[\(\{]\S*)','tokens');

for i = 1:length(props)
	% clean up quotes
	val = regexprep( props{i}(2),'["'']','' );
	% drop the % notation when given
	prop = regexprep(props{i}(1),'%','');
	% drop the '-' notation when given
	prop = regexprep(prop,'-','');
	
	if strcmp(lower(prop),'bus1')
		Phases=regexp(val,'\.','split');
		Phases=length([Phases{:}])-1;
		if Phases>0
			obj.Phases = Phases;
		end
	end
	
	numval = rpncal(val);
	if isempty(numval)
		obj.(prop) = val;
	else
		obj.(prop) = numval;
	end
end
end

function pathtofile = WriteDSS( dsscircuit, filename, splitFileFlag, savepath, commands)
% Write OpenDSS circuit to a single file or a set of files with each
% component stored in a file.
% Inputs:
%			dsscircuit: circuit object created from dssconversion function
%			filename: (optional) default: circuit's name. will be used as name for main circuit file and prefix for component files (e.g. [filename]_line.dss )
%			splitFileFlag: (optional) default: 0. Write data to multiple files with each component on a seperate file besides the main one.
%			savepath: (optional) relative/absolute path to save files. If folder doesn't exist, create one.
%			commands: additional commands
% Output:
%			pathtofile: path to main opendss file generated (useful for running OpenDSS Simulation in Matlab)

% Process inputs
if isstruct(dsscircuit)
	c = dsscircuit;
elseif ~isempty(strfind(class(dsscircuit),'dss'))
	splitFileFlag = 0;
	c.(class(dsscircuit)) = dsscircuit;
else
	error('Invalid data type for dsscircuit');
end
% remove wrong fields
cfields = fieldnames(dsscircuit);
for i = 1:numel(cfields)
	if isempty(strfind(class(dsscircuit.(cfields{i})), 'dss')) && ~strcmpi(cfields{i},'buslist') && ~strcmpi(cfields{i},'basevoltages')
		c= rmfield(c,cfields{i});
	end
end

if ~exist('splitFileFlag','var')
	splitFileFlag = 1;
end

headerfooterflag = 1;

if ~exist('savepath','var') || isempty(savepath)
	savepath = [pwd];
    if ~exist(savepath,'dir'), mkdir(savepath); end
else
	% handle relative path
	if(~strcmp(savepath(1:2),'\\') && savepath(2)~=':')
		savepath = [pwd '/' savepath];
	end
	
	% create folder if it doesn't exist
	if exist(savepath,'dir') < 1
		mkdir(savepath);
	end
end

% handle filename. Use circuit name if not specified
if exist('filename','var') && ~isempty(filename)
	if strfind(filename,'.dss'), 
		fname = filename(1:strfind(filename,'.dss')-1); 
	else
		fname = filename;
	end
else
	if ~isfield(c,'circuit') 
		warning('dsswrite:circuitUndefined','The input data doesn''t contain a circuit object! Check and make sure this is what you want.\nOpenDSS will not be able to load this file by itself.');
		fname = 'newdssfile';
		headerfooterflag = 0;
	else
		fname = c.circuit.Name;
	end
end

if exist([savepath '/' fname '.dss'],'file')
	try
		fnn = [savepath '/' fname '.dss'];
		fid = fopen(fnn,'w');
		fclose(fid);
	catch
		fname = [fname '_' datestr(now,'YYYYMMDDhhmmss')];
	end
end

% Writing out
s = '';
if isfield(c,'circuit')
	s = char(c.circuit);
end
fn = fieldnames(c);
% bus list isn't a class like the others
ind = strcmp(fn,'buslist')|strcmp(fn,'basevoltages');
if(any(ind))
	fn(ind) = [];
	if isfield(c,'buslist')
		buslist = c.buslist;
	end
end
% arrange them in the right order to print out (e.g. linecode should be
% defined before line)
classes = {'wiredata','linegeometry','linecode','line','loadshape','tshape','tcc_curve','reactor','fuse','transformer','regcontrol','capacitor','capcontrol','xycurve','pvsystem','InvControl','storage','storagecontroller','load','generator','swtcontrol','monitor','energymeter'};
classes(~ismember(classes,fn)) = [];

classes2 = {'switch','recloser'};
classes2(~ismember(classes2,fn)) = [];

fn = [classes, setdiff(fn,[classes classes2])', classes2];

% open main file
if isempty(strfind(fname,'.dss')), fname2 = [fname '.dss']; end
pathtofile = [savepath '/' fname2];
fidmain = fopen(pathtofile, 'w');
if fidmain < 1
	splitFileFlag = 0;
	fname2 = ['dss_' datestr(now,'YYYYMMDDhhmmss')];
	pathtofile = [savepath '/' fname2 '.dss'];
	fidmain = fopen(pathtofile, 'w');
end

if ~splitFileFlag
	% add each device class to circuit string
	for i = 1:length(fn)
		if strcmp('circuit',fn{i}), continue; end;
        if ~isempty(c.(fn{i}))
            s = [s char(c.(fn{i}))];
        end
	end
else
	% write files for all devices
	
	for i = 1:length(fn)
		
		%device filename
		dfn = [fname '_' fn{i} '.dss'];
		
		if strcmp('circuit',fn{i}), continue; end;
		% open file for writing
		fid = fopen([savepath '/' dfn], 'w');
		if(fid==-1), error('dsswrite:openfailed','Failed to open output file %s for writing!\nRemember to close open files before overwriting.',[savepath '/' dfn]); end
		s_ = char(c.(fn{i}));
		try
			fwrite(fid, s_);
			fclose(fid);
		catch err
			warning('dsswrite:openfiles','Remember to close files before overwriting them!');
			rethrow(err);
		end
			
		% update main file with "Redirect" command
		s = sprintf('%s\n%s',s,['Redirect ' dfn]); 
	end
end

if(~isfield(c,'basevoltages'))
	c.basevoltages = [115, 12.47, 4.16, 0.48, 0.12];
end

if headerfooterflag
	cvs = sprintf('\n\n! Let DSS estimate the voltage bases\n%s%s\n%s\n',...
		'Set voltagebases=',mat2str(c.basevoltages),...
		'Calcvoltagebases     ! This also establishes the bus list');
if isfield(c,'InvControl')
	iter=sprintf('\nSet maxcontroliter=5000');
	s = [s cvs iter];
else
	s = [s cvs];
end

end

if(exist('buslist','var'))
	s = sprintf('%s\nBuscoords %s_%s.csv\n', s, fname, 'buscoords');
	sbl = [buslist.id num2cell(buslist.coord)]';
	sbl = sprintf('%s, %g, %g\n',sbl{:});
	try
		fid = fopen([savepath '/' fname '_buscoords.csv'], 'w');
		fwrite(fid,sbl);
		fclose(fid);
	catch err
		warning('dsswrite:openfiles','Remember to close files before overwriting them!');
		rethrow(err);
	end
end

if headerfooterflag
	% header
	h = sprintf('%s\n\n','Clear');
	
	% footer
% 	f = sprintf('\n\n%s\n%s\n\n%s\n%s\n%s\n',...
% 				'set maxiterations=100',...
% 				'solve mode=snapshot',...
% 				'show voltages LL Nodes',...
% 				'show powers kva elements',...
% 				'show taps',...
%                 'export voltages',...
%                 'export seqvoltages',... 
%                 'export powers kva',...
%                 'export p_byphase',...
%                 'export seqpowers');		
	f = sprintf('\n\n%s\n%s\n\n%s\n%s\n%s\n',...
				'set maxiterations=1000');
else
	h = '';
	f = '';
end

% write main file
if ~exist('commands','var')
	fwrite(fidmain, [h s f]);
else
	fwrite(fidmain, sprintf('%s %s \n %s',h,s,commands));
end
fclose(fidmain);

end