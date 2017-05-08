%New FR test code
clear
Full=tic;

feeder='ieee13';
%% Todo
%%Capacitor
% Right now we keep all capacitors. By doing this we also keep both buses
% of the line the capacitor is connected to. In the future 1) We can remove
% capacitors which are not voltage sensative. 2) Remove the 2nd bus of the
% line? Move cap to shunt bus that controls voltage

%%Make sure phases are in order
%Take code from old code!

%%Shunt terms of distribution transformer
%Remove distribution transformer by replacing shunt terms with cap

%%length - maybe use topo detection!

%% start code
if strcmp(feeder,'ieee13')
tic
circuit=dssparse('C:\Users\Zactus\Documents\OpenDSS\OpenDSS\IEEETestCases\13Bus\IEEE13Nodeckt_zack.dss');
toc

[ data1 ] = dssSimulation( circuit,[],1,[],[],0)
p='C:\Users\Zactus\Documents\OpenDSS\OpenDSS\IEEETestCases\13Bus\IEEE13Nodeckt_zack.dss';
o = actxserver('OpendssEngine.dss');
dssText = o.Text; dssText.Command = 'Clear'; cDir = pwd;
dssText.Command = ['Compile "' p '"'];
dssCircuit = o.ActiveCircuit;
circuit.buslist_orig=circuit.buslist;
circuit.buslist.id=dssCircuit.AllBUSNames;
circuit.buslist.coord=zeros(length(circuit.buslist.id),2);
delete(o);
% circuit.linecode(end).cmatrix=[];
% circuit.linecode(end-1).cmatrix=[];
elseif strcmp(feeder,'Alpine')
circuit=load('c:\users\zactus\gridIntegration\results\Alpine_s1b_refined.mat');
circuit=circuit.c;

for j=1:length(circuit.pvsystem)
if isempty(regexp(circuit.pvsystem(j).bus1,'\.','match')) %single phase
circuit.pvsystem(j).kV=circuit.pvsystem(j).kV/sqrt(3);
end
end
% circuit=rmfield(circuit,'pvsystem');

[ data1 ] = dssSimulation( circuit,[],1,[],[],0)
end
% criticalBuses=circuit.buslist.id(round((length(circuit.buslist.id)-2)*rand(kk,1)+1));
% criticalBuses=unique(criticalBuses);
% while length(criticalBuses)<kk
% 	cnode_tmp=circuit.buslist.id(round((length(circuit.buslist.id)-2)*rand(kk,1)+1));
% 	criticalBuses=[criticalBuses;cnode_tmp];
% 	criticalBuses=unique(criticalBuses);
% end

criticalBuses=circuit.buslist.id([11,13,14]);

for ii=1:length(criticalBuses)
	criticalBuses(ii)=regexprep(criticalBuses(ii),'-','_');
end


%% check which variables we have and set flags
Flag=isfield(circuit,{'load','pvsystem','capacitor','transformer','regcontrol','capcontrol'});

%% Get Ybus and nameing vectors
if Flag(2)
	circuit_woPV=rmfield(circuit,'pvsystem');
else
	circuit_woPV=circuit;
end
if Flag(1);
	circuit_woPV=rmfield(circuit_woPV,'load');
end
circuit_woPV.line(end-1).Length=1;
circuit_woPV.linecode(end).Units='ft';
circuit_woPV.linecode(end).Cmatrix='[0]';
%load the circuit and generate the YBUS
p = WriteDSS(circuit_woPV,'test',0,'c:\users\zactus\FeederReductionRepo'); o = actxserver('OpendssEngine.dss');
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

busnames=regexp(dssCircuit.YNodeOrder,'\.','split');
YbusOrderVect=[busnames{:}]'; YbusOrderVect(find(cellfun('length',YbusOrderVect)==1))=[];
YbusPhaseVect=[busnames{:}]'; YbusPhaseVect(find(cellfun('length',YbusPhaseVect)>1))=[]; YbusPhaseVect=str2double(YbusPhaseVect);
Ycomb=strcat(YbusOrderVect,'.', num2str(YbusPhaseVect));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vahids voodoo to get Ybase
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
volt_base=volt_base(cell2mat(values(mapObj,valueSet)))

power_base=1000000;
Zbase=volt_base*volt_base'/power_base;
Ybase=(1./Zbase);

Ybus_pu=Ybus./Ybase;
%% %Here we generate the list of node numbers in the order of the Ybus
buslist=dssCircuit.AllBusNames;
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

%% get topo
topo=zeros(max(Node_number),4);
generation{1,1}=[];clear generation;  generation{1,1}=1; generation{1,4}=[];generation{1,5}=0;
parent=1;
topo(parent,1)=parent;
[topo,generation]=topology_detect_large(topo,generation,Ybus,parent,Node_number);
topo_view=topo;
topo_view(find(topo_view(:,1)==0)',:)=[];
c_new=0;

%% get transformers
if Flag(4)
	if Flag(5)
		rgc=circuit.regcontrol;
		for ii=1:length(rgc)
			trfName=rgc(ii).transformer;
			Trf_Ind(ii)=find(ismember({circuit.transformer{:}.Name},trfName));
% 			circuit.transformer(Trf_Ind).Name=['REG_' circuit.transformer(Trf_Ind).Name];
% 			circuit.regcontrol(ii).transformer=[circuit.transformer(Trf_Ind).Name];
			buses=circuit.transformer(Trf_Ind(ii)).buses;
			busNames=strtok(buses,'\.');
			criticalNumbers=[criticalNumbers; find(ismember(lower(buslist),lower(busNames)))];
		end
	end
	criticalNumbers=unique(criticalNumbers);
	Critical_buses=buslist(criticalNumbers);
end

% Delete all Xfrmrs not Reg or Substation
% % NonVRind=find(cellfun(@isempty,regexp(circuit.transformer(:).Name,'REG_','match')));
NonSubind=find(cellfun(@isempty,regexp(circuit.transformer(:).Name,'Sub','match')));
RemInd=find(~ismember(1:length(circuit.transformer),Trf_Ind));
RemTrans=RemInd(find(ismember(RemInd,NonSubind)));
circuit.transformer(RemTrans)=[];

%% get capacitors
%This section is to keep the capacitors in the grid by
%adding the nodes that they are connected to.
bus1=[];

if Flag(3)
	tic
	fprintf('\nAdding Capacitor Nodes: ')
	
	cap=circuit.capacitor;
	for ii=1:length(cap)
		bus1{ii}=cap(ii).bus1;
	end
	if ~iscolumn(bus1)
		bus1=bus1';
	end
	bus1=strtok(bus1,'\.'); CapNum=find(ismemberi(buslist,bus1));
	criticalNumbers=[criticalNumbers; CapNum];
	
	criticalBuses=buslist(criticalNumbers);
	clear bus1
	
	if Flag(6)
		capcon=circuit.capcontrol;
		for ii=1:length(capcon)
			buses=regexp(capcon(ii).Element,'\.','split');
			if strcmp(buses{1},'line')
				LineNo=find(ismember({circuit.line{:}.Name},buses{2}));
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
%% get load
%here we preserve by load type!
if Flag(1)
	fprintf('Initialize loads: ')
	tic
	ld=circuit.load;
	loadCurrents=zeros(length(YbusOrderVect),1);
	for j=1:length(ld)
		loadCurrent=(ld(j).Kw+1i*ld(j).Kvar);
		
		if isempty(regexp(ld(j).bus1,'\.','match')) %3 phase
			loadCurrents(find(strcmpi(YbusOrderVect,ld(j).bus1)==1),ld(j).model)=loadCurrent/3/(ld(j).kV/sqrt(3));
		elseif length(regexp(ld(j).bus1,'\.','match'))>1 %2 phase
			name=regexp(ld(j).bus1,'\.','split');
			numPhases=length(name)-1;
			for ii=2:length(name)
				loadCurrents(find(ismember(Ycomb,[name{1} '.' name{ii}])),ld(j).model)=loadCurrent/numPhases/(ld(j).kV/sqrt(3));
			end
		else %1 phase
			loadCurrents(find(ismember(Ycomb,ld(j).bus1)),ld(j).model)=loadCurrent./ld(j).kV;
		end
	end
	t_=toc;
	fprintf('%.2f sec\n',t_)
end

%% get PV
if Flag(2)
	fprintf('Initialize PV systems: ')
	tic
	pv=circuit.pvsystem;
	pvCurrents=zeros(length(YbusOrderVect),1);
	for j=1:length(pv)
		P=pv(j).pf*pv(j).kVA; Q=sqrt(pv(j).kVA^2-P^2);
		pvPower=(P+1i*Q);
		
		if isempty(regexp(pv(j).bus1,'\.','match')) %single phase
			pvCurrents(find(strcmpi(YbusOrderVect,pv(j).bus1)==1))=pvPower/3/(pv(j).kV/sqrt(3));
		elseif length(regexp(ld(j).bus1,'\.','match'))>1 %2 phase
			name=regexp(ld(j).bus1,'\.','split');
			for ii=2:length(name)
				pvCurrents(find(ismember(Ycomb,[name{1} '.' name{ii}])))=pvPower/length(2:length(name))/(pv(j).kV/sqrt(3));
			end
		else %3 phase
			pvCurrents(find(ismember(Ycomb,ld(j).bus1)))=pvPower/(pv(j).kV);
		end
	end
	t_=toc;
	fprintf('%.2f sec\n',t_)
end
%% First find the nodes the nodes that connect two or more critical nodes,
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
		%Here we find which parentso of CN(k) are also parents of CN(j),
		%return logic, find those indicices and return the parents of CN(k)
		VectorOfCommonParents=(generation{criticalNumbers(k),4}(find(ismember(cell2mat(generation(criticalNumbers(k),4)),cell2mat(generation(criticalNumbers(j),4))))));
		%find the node with the greatest distance (i.e. closest to each
		%node)
		[~,NewCriticalNumber]=max(cell2mat(generation(VectorOfCommonParents,5))); %Find the closest common point based on distance to substation (i.e. generation(_,5)
		New_criticalNumbers=[New_criticalNumbers, VectorOfCommonParents(NewCriticalNumber)];
	end
end
criticalNumbers=vertcat(criticalNumbers, New_criticalNumbers');

%Add sourcebus to CN
criticalNumbers=[criticalNumbers; find(ismember(lower(buslist),'sourcebus'))];

%Add all substation equipment to list of critical buses
for ii=1:length(buslist)
	if ~(isempty(cell2mat(regexp(buslist(ii),'sub'))))
		criticalNumbers=[criticalNumbers; ii];
	end
end

criticalNumbers=unique(criticalNumbers); %get rid of repeat connections
criticalBuses=buslist(criticalNumbers);

t_=toc;
fprintf('time elapsed %f\n',t_)

%% Actual Reduction Part
tic
fprintf('\nActual reduction: ')

Zbus=inv(Ybus);
BusOrder_new=find(ismember(lower(YbusOrderVect),criticalBuses));
Zbus_new=Zbus(BusOrder_new,BusOrder_new);
Ybus_reduced=inv(Zbus_new);
Ybus_reduced(find(abs(real(Ybus_reduced))<1E-9))=0;
Ybus_reduced(find(abs(imag(Ybus_reduced))<1E-9))=0;
buslist=criticalBuses;
YbusOrderVect=YbusOrderVect(BusOrder_new);
YbusPhaseVect=YbusPhaseVect(BusOrder_new);
Ycomb=Ycomb(BusOrder_new);

t_=toc;
fprintf('%.2f sec\n',t_)

% Load_init=sum(cell2mat(circuit.load(:).kW)+1i*cell2mat(circuit.load(:).kVar))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% new topo
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


%% get rid fo classes we no longer care aboot
names=fieldnames(circuit);
keepFields={'load','buslist','line','circuit','capcontrol','transformer','capacitor','basevoltages','regcontrol','pvsystem','reactor'};
names=names(find(~ismember(names,keepFields)));
for ii=1:length(names)
	circuit=rmfield(circuit,names{ii});
end

%% Rewrite lines
% Lengths(:,1)=strtok(circuit.line(:).bus1,'.');
% Lengths(:,2)=circuit.line(:).length;
% Lengths(:,3)=circuit.line(:).length;

trfBus=strtok(reshape([circuit.transformer{:}.buses],2,[]),'\.');

circuit.line=dssline;
count=0;
for ii=1:length(generation(:,1))
	Connected=cell2mat(generation(ii,2));
	Bus1=buslist(cell2mat(generation(ii,1)));
	Bus1Ind=find(ismember(lower(YbusOrderVect),Bus1));
	
	for jj=1:length(Connected)
		Bus2=buslist(Connected(jj));
		Bus2Ind=find(ismember(lower(YbusOrderVect),Bus2));
		
		%Make sure it is not VR
		if any(sum(ismember(lower(trfBus),[lower(Bus1); lower(Bus2)]))>1)
			continue
		else
			count=count+1;
			
			circuit.line(count)=dssline;
			circuit.line(count).Name=[char(Bus1) '_' char(Bus2)];
					circuit.line(count).Units='mi';
			circuit.line(count).R1=[];
			circuit.line(count).R0=[];
			circuit.line(count).X0=[];
			circuit.line(count).X1=[];
			circuit.line(count).C0=[];
			circuit.line(count).C1=[];
			
			Bus1IndMod=Bus1Ind(find(ismember(YbusPhaseVect(Bus1Ind),YbusPhaseVect(Bus2Ind))));
			circuit.line(count).bus1=[char(Bus1) '.' strjoin(arrayfun(@(x) num2str(x),YbusPhaseVect(Bus1IndMod),'UniformOutput',false),'.')];
			circuit.line(count).bus2=[char(Bus2) '.' strjoin(arrayfun(@(x) num2str(x),YbusPhaseVect(Bus2Ind),'UniformOutput',false),'.')];
			circuit.line(count).Phases=length(Bus2Ind);
			circuit.line(count).Length=1;
			
			ImpMat=-inv(Ybus_reduced(Bus1IndMod,Bus2Ind))./circuit.line(count).Length;
			
			if length(ImpMat)==1
				circuit.line(count).Rmatrix=['(' num2str(real(ImpMat(1,1))) ')'];
				circuit.line(count).Xmatrix=['(' num2str(imag(ImpMat(1,1))) ')'];
			elseif length(ImpMat)==2
				circuit.line(count).Rmatrix=['(' num2str(real(ImpMat(1,1))) '|' num2str(real(ImpMat(2,1:2)))  ')'];
				circuit.line(count).Xmatrix=['(' num2str(imag(ImpMat(1,1))) '|' num2str(imag(ImpMat(2,1:2)))  ')'];
			else
				circuit.line(count).Rmatrix=['(' num2str(real(ImpMat(1,1))) '|' num2str(real(ImpMat(2,1:2))) '|' num2str(real(ImpMat(3,1:3))) ')'];
				circuit.line(count).Xmatrix=['(' num2str(imag(ImpMat(1,1))) '|' num2str(imag(ImpMat(2,1:2))) '|' num2str(imag(ImpMat(3,1:3))) ')'];
			end
			
		end
	end
end

%update capcontrol to match updated line names
if isfield(circuit,'capcontrol')
	for ii=1:length(capConBuses1)
		Numms1=find(ismemberi(strtok(circuit.line(:).bus1,'.'),char(capConBuses1{ii})));
		Numms2=find(ismemberi(strtok(circuit.line(:).bus2,'.'),char(capConBuses2{ii})));
		lineInd=Numms1(find(ismemberi(Numms1,Numms2)));
		circuit.capcontrol(ii).Element=['line.' circuit.line(lineInd).Name];
	end
end
%% Rewrite loads
if Flag(1)
	fprintf('Re-writing PV systems: ')
	tic
	
	%Get buses downstream of Vreg to allocate Vbase to load and PV
	trfDown=trfBus(:,2);
	[~,trfDownInd]=ismember(lower(trfDown),lower(buslist));
	[trfDownInd, IA]=unique(trfDownInd); trfDown=trfDown(IA);
	trf_kV=cell2mat(reshape([circuit.transformer{:}.Kv],2,[])); trf_kV(1,:)=[];
	trf_phase=[circuit.transformer{:}.Phases];
	trf_phase(find(trf_phase==3))=4; trf_phase(find(trf_phase==1))=3; trf_phase(find(trf_phase==4))=1;
	trf_kV=trf_kV.*sqrt(trf_phase);  trf_kV=trf_kV(IA);
	
	%get load
	newLoadCurrents=Ybus_reduced*Zbus(BusOrder_new,:)*loadCurrents;
	
	count=0;
	circuit.load=dssload;
	for loadModel=1:size(newLoadCurrents,2)
		WriteLoads=find(newLoadCurrents(:,loadModel)>0.01);
		for ii=1:length(WriteLoads)
			count=count+1;
			circuit.load(count)=dssload;
			circuit.load(count).Name=['Load ' char(Ycomb(WriteLoads(ii))) '_model_' num2str(loadModel)];
			circuit.load(count).phases=1;
			
			busNum=find(ismemberi(buslist,YbusOrderVect(WriteLoads(ii))));
			MatchInd=find(ismemberi(trfDownInd,[cell2mat(generation(busNum,4)); busNum]));
			circuit.load(count).kV=trf_kV(MatchInd(end))/sqrt(3);
			
			circuit.load(count).bus1=Ycomb(WriteLoads(ii));
			circuit.load(count).Kw=real(newLoadCurrents(WriteLoads(ii),loadModel))*circuit.load(count).kV;
			circuit.load(count).KVAr=imag(newLoadCurrents(WriteLoads(ii),loadModel))*circuit.load(count).kV;
			circuit.load(count).model=loadModel;
			circuit.load(count).kVA=[];
			circuit.load(count).pf=[];
		end
	end
	t_=toc;
	fprintf('%.2f sec\n',t_)
end
% Load_end=sum(cell2mat(circuit.load(:).kW)+1i*cell2mat(circuit.load(:).kVar))
%% rewrite PV
if Flag(2)
	fprintf('Re-writing PV systems: ')
	tic
	
	%updated PV
	newPVCurrents=Ybus_reduced*Zbus(BusOrder_new,:)*pvCurrents;
	count=0;
	circuit.pvsystem=dsspvsystem;
	WritePVs=find(newPVCurrents(:)>0.01);
	for ii=1:length(WritePVs)
		count=count+1;
		circuit.pvsystem(count)=dsspvsystem;
		circuit.pvsystem(count).Name=['PV_' char(Ycomb(WritePVs(ii)))];
		circuit.pvsystem(count).phases=1;
		circuit.pvsystem(count).irradiance=1;
		circuit.pvsystem(count).Temperature=25;
		
		busNum=find(ismemberi(buslist,YbusOrderVect(WritePVs(ii))));
		MatchInd=find(ismemberi(trfDownInd,[cell2mat(generation(busNum,4)); busNum]));
		circuit.pvsystem(count).kV=trf_kV(MatchInd(end))/sqrt(3);
		
		circuit.pvsystem(count).bus1=Ycomb(WritePVs(ii));
		circuit.pvsystem(count).Pmpp=sqrt((real(newPVCurrents(WritePVs(ii)))*circuit.pvsystem(count).kV)^2+(imag(newPVCurrents(WritePVs(ii)))*circuit.pvsystem(count).kV)^2);
		circuit.pvsystem(count).kVA=sqrt((real(newPVCurrents(WritePVs(ii)))*circuit.pvsystem(count).kV)^2+(imag(newPVCurrents(WritePVs(ii)))*circuit.pvsystem(count).kV)^2);
		circuit.pvsystem(count).pf=(real(newPVCurrents(WritePVs(ii)))*circuit.pvsystem(count).kV)/circuit.pvsystem(count).kVA;
	end
	t_=toc;
	fprintf('%.2f sec\n',t_)
end

%% Update circuit info
circuit.circuit.Name=[circuit.circuit.Name '_Reduced'];
tmp=find(ismemberi(lower(circuit.buslist.id),buslist));
circuit.buslist.id=circuit.buslist.id(tmp);
circuit.buslist.coord=circuit.buslist.coord(tmp,:);

t_reduct=toc(Full);
if Flag(2)
	circuit_woPV=rmfield(circuit,'pvsystem');
else
	circuit_woPV=circuit;
end
if Flag(1);
	circuit_woPV=rmfield(circuit_woPV,'load');
end
p = WriteDSS(circuit_woPV,'Test',0,'c:\users\zactus\FeederReductionRepo\tmp');
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

busnames=regexp(dssCircuit.YNodeOrder,'\.','split');
YbusOrderVect2=[busnames{:}]'; YbusOrderVect2(find(cellfun('length',YbusOrderVect2)==1))=[];
YbusPhaseVect2=[busnames{:}]'; YbusPhaseVect2(find(cellfun('length',YbusPhaseVect2)>1))=[]; YbusPhaseVect2=str2double(YbusPhaseVect2);
Ycomb2=strcat(YbusOrderVect2,'.', num2str(YbusPhaseVect2));

[ data ] = dssSimulation( circuit,[],1,[],[],0)
Keep=find(ismember(data1.nodeName,data.nodeName));
figure;plot(1:length(Keep),data.Voltage,'r*',1:length(Keep),data1.Voltage(Keep),'b*')

Order=[]; Prev=[];
for ii=1:length(YbusOrderVect)
Current=find(ismember(YbusOrderVect2,YbusOrderVect(ii)));
if ismember(Prev,Current)
Prev=Current;
continue
end
Order=[Order; Current];
Prev=Current;
end
YbusOrderVect_New=YbusOrderVect2(Order);
Ybus_reduced_ordered=Ybus_reduced(Order,Order);
Ydiff=Ybus_regenerated-Ybus_reduced_ordered;
[maxA,ind] = max(abs(real(Ydiff(:))))
[maxA,ind] = max(abs(imag(Ydiff(:))))
[m,n] = ind2sub(size(Ydiff),ind)
YbusOrderVect([m n])
