%New FR test code
clear
Full=tic;
close all

feeder='ieee13';
%% Todo
%%Capacitor
% Right now we keep all capacitors. By doing this we also keep both buses
% of the line the capacitor is connected to. In the future 1) We can remove
% capacitors which are not voltage sensative. 2) Remove the 2nd bus of the
% line? Move cap to shunt bus that controls voltage
%Also add in shunt terms
%Also, don't need CAP element (just need cap buses) MAke element the cap
%bus and any connected line

%%switches
%When switch == yes, the Ybus is wrong.
%When switch == no, the length is wrong






%% Notes on validation
% For NO reduction, just rewrite
% 1) Ybus is the same (1E-6 error) before and after reduction
% 2) Total kVAr differs by same b4 and after reduction


%% start code
if strcmp(feeder,'ieee13')
	tic
% 		p='C:\Users\Zactus\FeederReductionRepo\8500-Node\Master.dss';
		p='C:\Users\Zactus\Documents\OpenDSS\OpenDSS\IEEETestCases\13Bus\IEEE13Nodeckt_zack.dss';
% 		p='C:\Users\Zactus\Documents\OpenDSS\OpenDSS\IEEETestCases\123Bus\IEEE123Master.dss';
% 		p='C:\Users\Zactus\Documents\OpenDSS\OpenDSS\IEEETestCases\34Bus\ieee34_1.dss';
% 	circuit=dssparse('C:\Users\Zactus\Documents\OpenDSS\OpenDSS\IEEETestCases\13Bus\IEEE13Nodeckt_zack.dss');
% 	if ~exist('c:\users\zactus\FeederReductionRepo\ieee8500.mat')
		circuit=dssparse(p);
% 		save('c:\users\zactus\FeederReductionRepo\ieee8500.mat','circuit')
% 	else
% 		load('c:\users\zactus\FeederReductionRepo\ieee8500.mat')
% 	end
	toc
	
% 	orig=tic;
% 	[ data1 ] = dssSimulation( circuit,[],1,[],[],0);
% 	FullTime=toc(orig);
	
	o = actxserver('OpendssEngine.dss');
	dssText = o.Text; dssText.Command = 'Clear'; cDir = pwd;
	dssText.Command = ['Compile "' p '"'];
	dssCircuit = o.ActiveCircuit;
	circuit.buslist.id=dssCircuit.AllBUSNames;
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
criticalBuses=circuit.buslist.id([1:1000]);
% criticalBuses=circuit.buslist.id;


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


for ii=1:length(criticalBuses)
	criticalBuses(ii)=regexprep(criticalBuses(ii),'-','_');
end


%% check which variables we have and set flags
[Flag]=isfield(circuit,{'load','pvsystem','capacitor','transformer','regcontrol','capcontrol','reactor'});
circuit_orig=circuit;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Y-BUS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



%load the circuit and generate the YBUS
p = WriteDSS(circuitBase,'test',0,'c:\users\zactus\FeederReductionRepo'); o = actxserver('OpendssEngine.dss');

% p = WriteDSS(circuit,'test',0,'c:\users\zactus\FeederReductionRepo'); o = actxserver('OpendssEngine.dss');
dssText = o.Text; dssText.Command = 'Clear'; cDir = pwd;
dssText.Command = ['Compile "' p '"']; dssCircuit = o.ActiveCircuit;

dssSolution = dssCircuit.Solution;
dssSolution.MaxControlIterations=100;
dssSolution.MaxIterations=100;
dssSolution.InitSnap; % Initialize Snapshot solution
dssSolution.dblHour = 0.0;
dssSolution.Solve;

Ybus=dssCircuit.SystemY;
ineven=2:2:length(Ybus); inodd=1:2:length(Ybus);
Ybus=Ybus(inodd)+1i*Ybus(ineven); Ybus=reshape(Ybus,sqrt(length(Ybus)),sqrt(length(Ybus)));

[YbusOrderVect, YbusPhaseVect]=strtok(dssCircuit.YNodeOrder,'\.');
YbusPhaseVect=str2num(cell2mat(strrep(YbusPhaseVect,'.','')));
Ycomb=dssCircuit.YNodeOrder;

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
keySet = lower(Ycomb);
valueSet = lower(dssCircuit.AllnodeNames);
mapObj = containers.Map(keySet,1:length(keySet));
volt_base=volt_base(cell2mat(values(mapObj,valueSet)));

power_base=1;
Zbase=volt_base*volt_base'/power_base;
Ybase=(1./Zbase);

circuit_orig.Ybus=Ybus;
% Ybus=Ybus./Ybase;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Buslist    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

buslist=dssCircuit.AllBusNames;
%get rid of this using key vals!!!!!!!!!!!!!!!
for ii=1:length(buslist)
	Ind=find(strcmpi(buslist(ii),YbusOrderVect))';
	Node_number(Ind)=ii;
end
for j=1:length(criticalBuses)
	criticalNumbers(j)=find(strcmpi(buslist,lower(criticalBuses(j))));
end
criticalNumbers=unique(criticalNumbers);
criticalBuses=buslist(criticalNumbers);

if ~iscolumn(criticalNumbers)
	criticalNumbers=criticalNumbers';
end

%Add sourcebus


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     TOPO    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topo=zeros(max(Node_number),4);
generation{1,1}=[];clear generation;  generation{1,1}=1; generation{1,4}=[];generation{1,5}=0;
parent=1;
topo(parent,1)=parent;
[topo,generation]=topology_detect_large(topo,generation,Ybus,parent,Node_number);
topo_view=topo;
topo_view(find(topo_view(:,1)==0)',:)=[];
c_new=0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    X-frmrs    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Flag(4)
	if Flag(5)
		rgc=circuit.regcontrol;
		for ii=1:length(rgc)
			trfName=rgc(ii).transformer;
			Trf_Ind(ii)=find(ismemberi({circuit.transformer{:}.Name},trfName));
			buses=circuit.transformer(Trf_Ind(ii)).buses;
			busNames=strtok(buses,'\.');
			criticalNumbers=[criticalNumbers; find(ismemberi(buslist,busNames))];
		end
	end
	criticalNumbers=unique(criticalNumbers);
	criticalBuses=buslist(criticalNumbers);
end

% Delete all Xfrmrs not Reg or Substation
% % NonVRind=find(cellfun(@isempty,regexp(circuit.transformer(:).Name,'REG_','match')));
NonSubind=find(cellfun(@isempty,regexp(lower(circuit.transformer(:).Name),'sub','match')));
RemInd=find(~ismemberi(1:length(circuit.transformer),Trf_Ind));
RemTrans=RemInd(find(ismemberi(RemInd,NonSubind)));
circuit.transformer(RemTrans)=[];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Capacitors    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This section is to keep the capacitors in the grid by adding the nodes that they are connected to.

bus1=[];

if Flag(3)
	tic
	fprintf('\nAdding Capacitor Nodes: ')
	
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
	
	t_=toc;
	fprintf('time elapsed %f\n',t_)
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     LOADS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%here we preserve by load type!


if Flag(1)
	fprintf('Initialize loads: ')
	tic
	ld=circuit.load;
	S_LD=zeros(length(YbusOrderVect),8);
% 	kV_LD=zeros(length(YbusOrderVect),length(YbusOrderVect));
	Ra=(1/sqrt(3))*(cosd(-30)+1i*sind(-30));
	Rb=(1/sqrt(3))*(cosd(30)+1i*sind(30));
	
	for j=1:length(ld)
		s_ld=(ld(j).Kw+1i*ld(j).Kvar);
		
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
% 			kV_LD(Ind,Ind)=ld(j).kV;
		end
		
	end
% 	kV_LD(find(YbusPhaseVect==2))=kV_LD(find(YbusPhaseVect==2)).*exp(1i*120*(pi/180));
% 	kV_LD(find(YbusPhaseVect==3))=kV_LD(find(YbusPhaseVect==3)).*exp(-1i*120*(pi/180));
	
	t_=toc;
	fprintf('%.2f sec\n',t_)
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     PV    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Flag(2)
	fprintf('Initialize PV systems: ')
	tic
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
	t_=toc;
	fprintf('%.2f sec\n',t_)
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     TOPOGRAPHY NODES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this will now also become a critical node. This is done by looking at the
%grandparents of each critical node and seeing if there are common nodes
%between two critical nodes. The closest common node is the one that is
%kept
tic
fprintf('\nGetting topogrophical critical Nodes (connection): ')
nn=length(criticalNumbers);
New_criticalNumbers=[];
for k=1:nn
	for j=k+1:nn
		if ~isempty(generation{criticalNumbers(k),4})
		%Here we find which parentso of CN(k) are also parents of CN(j),
		%return logic, find those indicices and return the parents of CN(k)
		VectorOfCommonParents=(generation{criticalNumbers(k),4}(find(ismemberi(cell2mat(generation(criticalNumbers(k),4)),cell2mat(generation(criticalNumbers(j),4))))));
		%find the node with the greatest distance (i.e. closest to each
		%node)
		[~,NewCriticalNumber]=max(cell2mat(generation(VectorOfCommonParents,5))); %Find the closest common point based on distance to substation (i.e. generation(_,5)
		% 		[buslist(criticalNumbers(k)) buslist(criticalNumbers(j)) buslist(VectorOfCommonParents(NewCriticalNumber))]
		New_criticalNumbers=[New_criticalNumbers, VectorOfCommonParents(NewCriticalNumber)];
		else
			break
		end
	end
end
criticalNumbers=vertcat(criticalNumbers, New_criticalNumbers');

%Add sourcebus to CN
criticalNumbers=[criticalNumbers; find(ismemberi(buslist,'sourcebus'))];

%Add all substation equipment to list of critical buses
criticalNumbers=[criticalNumbers; find(not(cellfun('isempty', strfind(lower(buslist), 'sub'))))];

% for ii=1:length(buslist)
% 	if ~(isempty(cell2mat(regexp(buslist(ii),'sub'))))
% 		criticalNumbers=[criticalNumbers; ii];
% 	end
% end
% Index = find(not(cellfun('isempty', strfind(lower(buslist), 'sub'))));

criticalNumbers=unique(criticalNumbers); %get rid of repeat connections
criticalBuses=buslist(criticalNumbers);

t_=toc;
fprintf('time elapsed %f\n',t_)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reduction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Actual Reduction Part
tic
fprintf('\nActual reduction: ')


Zbus=inv(Ybus);
YbusOrder_reduced=find(ismemberi(YbusOrderVect,criticalBuses));
Zbus_new=Zbus(YbusOrder_reduced,YbusOrder_reduced);
Ybus_reduced=inv(Zbus_new);
Ybus_reduced(find(abs(Ybus_reduced)<1E-8))=0;

%keep track of order
buslist=criticalBuses; YbusOrderVect=YbusOrderVect(YbusOrder_reduced); YbusPhaseVect=YbusPhaseVect(YbusOrder_reduced); Ycomb=Ycomb(YbusOrder_reduced);

W=Zbus_new\Zbus(YbusOrder_reduced,:);
W(find(abs(W)<1E-6))=0;

V2=diag(volt_base./sqrt(3));
	
if Flag(1)
S_LD_new=V2(YbusOrder_reduced,YbusOrder_reduced)*conj(W)*diag(1./diag(V2))*S_LD;
fprintf('\nDiff between load before and after reduction: %f kW\n',sum(sum(S_LD))-sum(sum(S_LD_new)))
end
if Flag(2)
S_PV_new=V2(YbusOrder_reduced,YbusOrder_reduced)*conj(W)*diag(1./diag(V2))*S_PV;
fprintf('\nDiff between pv before and after reduction: %f kW\n',sum(sum(S_PV))-sum(sum(S_PV_new)))
end

t_=toc;
fprintf('%.2f sec\n',t_)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REWRITE CIRCUIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   CLEAN-UP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names=fieldnames(circuit);
keepFields={'load','buslist','line','circuit','capcontrol','transformer','capacitor','basevoltages','regcontrol','pvsystem','reactor'};
names=names(find(~ismemberi(names,keepFields)));
for ii=1:length(names)
	circuit=rmfield(circuit,names{ii});
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   LINES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trfBus=strtok(reshape([circuit.transformer{:}.buses],2,[]),'\.');

Lines(:,1)=strtok(circuit.line(:).bus1,'.');
Lines(:,2)=strtok(circuit.line(:).bus2,'.');
Lines(:,3)=circuit.line(:).length;
Lines(:,4)=circuit.line(:).units;
[Lines(:,3),Lines(:,4)]=Convert_to_kft(Lines(:,4),Lines(:,3));

%find lines that are switches
Inds=[find(ismemberi(circuit.line(:).switch,'true')) find(ismemberi(circuit.line(:).switch,'yes'))];
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
		
		
		%Make sure it is not VR
		if any(sum(ismemberi(trfBus,[Bus1; Bus2]))>1)
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     CAPACITORS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get buses downstream of Vreg to allocate Vbase to load and PV
trfDown=trfBus(2,:);
[~,trfDownInd]=ismemberi(trfDown,buslist);
[trfDownInd, IA]=unique(trfDownInd); trfDown=trfDown(IA);
trf_kV=cell2mat(reshape([circuit.transformer{:}.Kv],2,[])); trf_kV(1,:)=[];
trf_phase=[circuit.transformer{:}.Phases];
trf_phase(find(trf_phase==3))=4; trf_phase(find(trf_phase==1))=3; trf_phase(find(trf_phase==4))=1;
trf_kV=trf_kV(IA); %trf_kV=trf_kV.*sqrt(trf_phase);

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
ShuntRows=find(abs(sum(Ybus_reduced,2))>1E-6);
CapRows=ShuntRows(find(~ismemberi(YbusOrderVect(ShuntRows),unique(trfBus))));
trfRows=[];capRows=[];
if Flag(4)
trfRows=[circuit.transformer{:}.buses];
end
if Flag(3)
capRows=[circuit.capacitor(:).bus1];
end
kill=[trfRows'; capRows; 'sourcebus'];
killRows=unique([find(ismemberi(Ycomb,kill)); find(ismemberi(YbusOrderVect,kill))]);
keepRows=CapRows(find(~ismemberi(CapRows,killRows)));

% %Make additional capacitors here to account for removed impedence from DT
for ii=1:length(keepRows)
	if ~Flag(7) && ii==1
	circuit.reactor(1)=dssreactor;
	else
	circuit.reactor(end+1)=dssreactor;
	end
	circuit.reactor(end).name=['addedShunt_' char(regexprep(Ycomb{keepRows(ii)},'\.','_'))];
	circuit.reactor(end).phases=1;
	circuit.reactor(end).bus1=char(Ycomb(keepRows(ii)));
	r_jx=1/sum(Ybus_reduced(keepRows(ii),:),2);
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

p = WriteDSS(circuitBase,'Test',0,'c:\users\zactus\FeederReductionRepo\tmp');
% p = WriteDSS(circuit,'Test',0,'c:\users\zactus\FeederReductionRepo\tmp');
fprintf('total reduction: %.2f sec\n',t_reduct)
o = actxserver('OpendssEngine.dss');
dssText = o.Text; dssText.Command = 'Clear'; cDir = pwd;
dssText.Command = ['Compile "' p '"']; dssCircuit = o.ActiveCircuit;

dssSolution = dssCircuit.Solution;
dssSolution.MaxControlIterations=100;
dssSolution.MaxIterations=100;
dssSolution.InitSnap; % Initialize Snapshot solution
dssSolution.dblHour = 0.0;
dssSolution.Solve;

Ybus_regenerated=dssCircuit.SystemY;
%Convert the Ybus_regenerated to a matrix
ineven=2:2:length(Ybus_regenerated); inodd=1:2:length(Ybus_regenerated);
Ybus_regenerated=Ybus_regenerated(inodd)+1i*Ybus_regenerated(ineven); Ybus_regenerated=reshape(Ybus_regenerated,sqrt(length(Ybus_regenerated)),sqrt(length(Ybus_regenerated)));

[YbusOrderVect2, YbusPhaseVect2]=strtok(dssCircuit.YNodeOrder,'\.');
YbusPhaseVect2=str2num(cell2mat(strrep(YbusPhaseVect2,'.','')));
Ycomb2=dssCircuit.YNodeOrder;

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

figure;plot(data.Dist,data.Voltage,'r*',data1.Dist(Keep),data1.Voltage(Keep),'b*')
legend('Reduced','Original')
title('Matching Buses','fontsize',14)
xlabel('Distance to substation [km]','fontsize',14)
ylabel('Bus Voltage [V pu]','fontsize',14)
figure;plot(data.Dist,data.Voltage,'r*',data1.Dist,data1.Voltage,'b*')
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
