function [outputdss] = FeederReduction_dss(pathToDSSFile,critical_nodes)

%Created by Zachary K. Pecenak on 12/14/2015

%The purpose of this function is to reduce an OpenDSS feeder to a few
%selected node. A user might want to reduce a feeder for several reasons,
%but mainly to reduce the system complexity for increased computational
%efficiency. This function works by taking in the openDSS feeder model and 
%reducing the impedences, loads, and PV to a few specified nodes. Note that
%nodes connecting two or more critical nodes are inherently also critical
%nodes.
%The critical nodes are the ones that the user wants to keep in the reduced
%circuit. Thdey are input as a cell of the node names in the circuit file
%or the node number 
%The circuit is the c file for the specific feeder setup, and can be found
%in the saved feeder setup file.
%weights are stored in circuit.weights

%Example input
% load('c:\users\zactus\gridIntegration\results\ValleyCenter_wpv_existing.mat')
% critical_nodes={'03551325','03553322','035560','03555704','0355948','03552613','03552637','03551382','03554401','03552641A'};
% or critical_nodes=c.buslist.id; this is the full feeder
% [c] = FeederReduction(critical_nodes,c);

% 
% %Check both inputs are met
% if nargin<1
% 	error('This requires two inputs, the desired nodes to keep and the feeder circuit')
% end

circuit=dssparse(pathToDSSFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% Initial cleanup of data

%Apparently some of the phases on the transfromers were wrong, so we are
%cleanign that up and correctign them to what they are supposed to be
if isfield(circuit,'transformer')
for ii=1:length(circuit.transformer)
	Phases=sum(ismember(char(circuit.transformer(ii).buses(2)),'\.'));
	if Phases==0
		circuit.transformer(ii).Phases=3;
	else
	circuit.transformer(ii).Phases=Phases;
	end
end
end
% % %Get linelengths from orig circuit
% % for ii=1:length(circuit.line)
% % 	bus1=regexp(circuit.line(ii).bus1,'\.','split');
% % 	bus1(find(cellfun('length',bus1')==1))=[];
% % 	bus2=regexp(circuit.line(ii).bus2,'\.','split');
% % 	bus2(find(cellfun('length',bus2')==1))=[];
% % 	LengthVect(ii,1)=bus1;
% % 	LengthVect(ii,2)=bus2;
% % 	LengthVect(ii,3)={circuit.line(ii).Length};
% % end

%make sure nodes are column
if ~iscolumn(critical_nodes)
	critical_nodes=critical_nodes';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% get Ybus from OpenDSS for feeder.
fprintf('\nGetting YBUS and buslist from OpenDSS: ')
tic
%remove load and PV from Ybus
circuit_woPV=rmfield(circuit,'pvsystem');
circuit_woPV=rmfield(circuit_woPV,'load');
%load the circuit and generate the YBUS
o = actxserver('OpendssEngine.dss');
dssText = o.Text; dssText.Command = 'Clear'; cDir = pwd;
dssText.Command = ['Compile "' pathToDSSFile '"']; dssCircuit = o.ActiveCircuit;
Ybus=dssCircuit.SystemY;

%Convert the Ybus to a matrix
ineven=2:2:length(Ybus); inodd=1:2:length(Ybus);  
Ybus=Ybus(inodd)+1i*Ybus(ineven); Ybus=reshape(Ybus,sqrt(length(Ybus)),sqrt(length(Ybus)));
Ybus=sparse(Ybus);
%get buslist in order of Ybus and rearrange
busnames=regexp(dssCircuit.YNodeOrder,'\.','split');
YbusOrderVect=[busnames{:}]'; YbusOrderVect(find(cellfun('length',YbusOrderVect)==1))=[];
YbusPhaseVect=[busnames{:}]'; YbusPhaseVect(find(cellfun('length',YbusPhaseVect)>1))=[]; YbusPhaseVect=str2double(YbusPhaseVect);
% %Here we generate the list of node numbers in the order of the Ybus
buslist=dssCircuit.AllBusNames;
Origbuslist=buslist;
 for ii=1:length(buslist)
Ind=find(strcmpi(buslist(ii),YbusOrderVect))';
Node_number(Ind)=ii;
 end

clear inodd ineven 

t_=toc;
fprintf('time elapsed %f',t_)

%Check to see that the critical nodes are in the circuit
if ischar(critical_nodes)
	critical_nodes={critical_nodes};
end
critical_nodes(find(~ismember(lower(critical_nodes),lower(buslist))))=[];
if ~any(ismember(lower(critical_nodes),lower(buslist)))
 	error('The following nodes arent in the circuit: \n%s', critical_nodes{~ismember(critical_nodes,buslist)})
end
delete(o)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% Get topogrophy %% written by Vahid R. Disfani
%Uncommented
tic
fprintf('\nGetting topogrophy: ')
topo=zeros(max(Node_number),4);
generation{1,1}=[];clear generation;  generation{1,1}=1; generation{1,4}=[];generation{1,5}=0;
parent=1;
topo(parent,1)=parent;
[topo,generation]=topology_detect_large(topo,generation,Ybus,parent,Node_number);
topo_view=topo;
topo_view(find(topo_view(:,1)==0)',:)=[];
c_new=0;

% % % treeplot(topo(:,2)')
% % % [x,y] = treelayout(topo(:,2)');
% % % for ii=1:length(x)
% % %     text(x(ii),y(ii),num2str(ii))
% % % end

t_=toc;
fprintf('time elapsed %f',t_)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% Get line lengths and buses
for ii=1:length(circuit.line)
	Lines{ii,1}=circuit.line(ii).Length;
	Lines{ii,2}=regexprep(circuit.line(ii).bus1,'\..*','');
	Lines{ii,3}=regexprep(circuit.line(ii).bus2,'\..*','');
end

%% get pv data and match with nodes
%Format the various PV things

fprintf('\nGetting PV and load of circuit: ')

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% get load data and match with nodes
%All of this is basically to just make the load 3 phase
%There is undoubtably a more efficient way
load=circuit.load;
LOAD=zeros(length(YbusOrderVect),1);
for ii=1:length(YbusOrderVect)
	LoadBusOrderVect{ii}=[char(YbusOrderVect{ii}) '.' num2str(YbusPhaseVect(ii))];
end

for j=1:length(load)
	if isempty(regexp(load(j).bus1,'\.','match'))
		LOAD(find(strcmpi(YbusOrderVect,load(j).bus1)==1))=load(j).Kw/3+i*load(j).Kvar/3;
	end
	LOAD(find(strcmpi(LoadBusOrderVect,load(j).bus1)==1))=load(j).Kw+i*load(j).Kvar;	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% match critical nodes with node numbers
for j=1:length(critical_nodes)
Critical_numbers(j)=find(strcmpi(buslist,critical_nodes(j)));
end
Orig_nodes=critical_nodes;
t_=toc;
fprintf('time elapsed %f',t_)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% Add capacitors + transformers to Critical Nodes
%This section is to keep the capacitors and the transformers in the grid by
%adding the nodes that they are connected to.
tic
fprintf('\nAdding Capacitor and Transformer Nodes: ')

%Capacitor
if isfield(circuit,'capacitor')
	cap=circuit.capacitor;
	for ii=1:length(cap)
		bus1{ii}=cap(ii).bus1;
	end
	if ~iscolumn(bus1)
		bus1=bus1';
	end
	bus1=regexp(bus1,'\.','split');
	bus1=[bus1{:}]'; bus1(find(cellfun('length',bus1)==1))=[];
	New_Node_Num=find(ismember(buslist,bus1));
	critical_nodes=[critical_nodes; bus1];
	Critical_numbers=[Critical_numbers'; New_Node_Num];
end
clear bus1
% fprintf('Number of Critical Nodes)
%capactior controller
if isfield(circuit,'capcontrol')
	capcon=circuit.capcontrol;
	for ii=1:length(capcon)
		buses=regexp(capcon(ii).Element,'\.','split');
		buses=regexp(buses{2},'\_','split');
		bus1{ii}=buses{1};
		bus2{ii}=buses{2};
	end
	bus1=regexp(bus1,'\.','split');
	bus1=[bus1{:}]'; bus1(find(cellfun('length',bus1)==1))=[];
	capConBuses1=bus1';
	bus2=regexp(bus2,'\.','split');
	bus2=[bus2{:}]'; bus2(find(cellfun('length',bus2)==1))=[];
	capConBuses2=bus2';
	New_Node_Num=reshape([find(ismember(lower(buslist),lower(bus1))), find(ismember(lower(buslist),lower(bus2)))],[],1);
	critical_nodes=[critical_nodes; bus1; bus2];
	Critical_numbers=[Critical_numbers; New_Node_Num];
end
clear bus1 bus2
%transformer
if isfield(circuit,'transformer')
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
	New_Node_Num=reshape([find(ismember(lower(buslist),lower(bus1))), find(ismember(lower(buslist),lower(bus2)))],[],1);
	critical_nodes=[critical_nodes; bus1; bus2];
	Critical_numbers=[Critical_numbers; New_Node_Num];
end

clear bus1 bus2 buses

t_=toc;
fprintf('time elapsed %f',t_)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% First find the nodes the nodes that connect two or more critical nodes,
%this will now also become a critical node. This is done by looking at the
%grandparents of each critical node and seeing if there are common nodes
%between two critical nodes. The closest common node is the one that is
%kept
tic
fprintf('\nGetting topogrophical critical Nodes (connection): ')
Critical_numbers=unique(Critical_numbers); %get rid of repeat connections
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
Critical_numbers=[Critical_numbers; 1];
Critical_numbers=unique(Critical_numbers); %get rid of repeat connections
critical_nodes=buslist(Critical_numbers);
% % % % treeplot(topo(:,2)')
% % % % [x,y] = treelayout(topo(:,2)');
% % % % for ii=1:length(x)
% % % %     text(x(ii),y(ii),num2str(ii))
% % % % end
t_=toc;
fprintf('time elapsed %f',t_)
fprintf('\n\nTotal critical nodes= %d\n',length(Critical_numbers))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% Anything extra
%Get line names to store lengths
lineBus1=regexp(circuit.line{ii}.bus1,'\.','split','once');	
lineBus2=regexp(circuit.line{ii}.bus2,'\.','split','once');
for ii=1:length(circuit.line)
	LineBus1{ii}=lineBus1{ii,1}(1);
	LineBus2{ii}=lineBus2{ii,1}(1);
end

%Store Variables that are usefull later
AllCN=critical_nodes;
OriginalYbusOrderVect=YbusOrderVect;
OriginalYbusPhaseVect=YbusPhaseVect;
Originalbuslist=buslist;


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
		PvMat=zeros(3,1);
		LoadMat=zeros(3,1);
		
		%Loop through buses that need to leave town
		for jj=1:length(EndBuses)
			
			% Find the Node indices corresponding to the end bus
			EndNodeInd=find(ismember(lower(YbusOrderVect),buslist(EndBuses(jj))));												% Find Indices of the individual end node
			
			% Find the indices of the phases that correspond to the end node and the
			% critical node
			EndPhasesInd=find(ismember(YbusPhaseVect(EndNodeInd),YbusPhaseVect(Ind)));												% Find the phases of that line, to keep only proper phases
			
			% Find the non-zero pv nodes - means there is PV on that node
			PvPhases=find(PV(EndNodeInd));
			
			% Find the phases of the CN that are also phases of the PV
			MatchPVandCN=find(ismember(YbusPhaseVect(Ind),YbusPhaseVect(EndNodeInd(PvPhases))));
			
			%Update PV and Load
			PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+PV(EndNodeInd(EndPhasesInd));			% Update Pv
			LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+LOAD(EndNodeInd(EndPhasesInd));	% update load
			
% 			if find(ismember('035543',buslist(EndBuses(jj))))
% 				stop=1;
% 			end
			
			%Update weights - Add the end node rows with PV to the CN nodes
			%which match the PV node phases
			Weights(Ind(MatchPVandCN),:)=Weights(Ind(MatchPVandCN),:)+Weights(EndNodeInd(PvPhases),:);																		% Update weights
		end
		
		LOAD(Ind)=LOAD(Ind)+LoadMat(YbusPhaseVect(Ind));
		PV(Ind)=PV(Ind)+PvMat(YbusPhaseVect(Ind));
		
		%store children to delete at end. If you delte in real tiem it will
		%mess with the numbering.
		Store=[Store; EndBuses];
		
	end
end
t_=toc;
fprintf('time elapsed %f',t_)

%Update store to reflect Ybus order
Store=unique(Store);
S=find(ismember(lower(YbusOrderVect),buslist(Store)));

%Reduce everything to collapse chitlins
Ybus(S,:)=[]; Ybus(:,S)=[];
Lines(find(ismember(lower(Lines(:,2)),buslist(Store))),:)=[];
Lines(find(ismember(lower(Lines(:,3)),buslist(Store))),:)=[];
buslist(Store)=[];
Reduction(1)=length(Store);
PV(S)=[];
LOAD(S,:)=[];
YbusOrderVect(S)=[];
YbusPhaseVect(S)=[];
Node_number=[];
for ii=1:length(buslist)
Ind=find(strcmpi(buslist(ii),YbusOrderVect))';
Node_number(Ind)=ii;
end
Weights(S,:)=[];


%Need to update critical numbers to busnames
Critical_numbers=find(ismember(buslist,critical_nodes));
Store=[];

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
	CNchildren=Critical_numbers(find(ismember(Critical_numbers,generation{Critical_numbers(j),3})));
	GenToRem=CNchildren;
	
	%remove branches that are attached to critical nodes that don't have
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
	PvMat=zeros(3,1);
	LoadMat=zeros(3,1);

	%Keep only corect phases
	for jj=1:length(AllGen)
		
		%Find the nodes that correspond to the bus being removed
		EndNodeInd=find(ismember(lower(YbusOrderVect),lower(buslist(AllGen(jj)))));												% Find Indices of the individual end node
		
		%find the indices end node phases that are also CN phases
		EndPhasesInd=find(ismember(YbusPhaseVect(EndNodeInd),YbusPhaseVect(Ind)));
		
		%Find the end nodes which have PV on them
		PvPhases=find(PV(EndNodeInd));
		
		%Match the CN Phases to update with the phases on end node that
		%have PV
		MatchPVandCN=find(ismember(YbusPhaseVect(Ind),YbusPhaseVect(EndNodeInd(PvPhases))));
		
		
		PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+PV(EndNodeInd(EndPhasesInd));			% Update Pv
		LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+LOAD(EndNodeInd(EndPhasesInd));	% update load
		
		%Update weights
		Weights(Ind(MatchPVandCN),:)=Weights(Ind(MatchPVandCN),:)+Weights(EndNodeInd(PvPhases),:);																	% Update weights
	end
	
	%Update PV and Load
	LOAD(Ind)=LOAD(Ind)+LoadMat(YbusPhaseVect(Ind));
	PV(Ind)=PV(Ind)+PvMat(YbusPhaseVect(Ind));
	
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
			PvMat=zeros(3,1);
			LoadMat=zeros(3,1);
			
			
			%Keep only corect phases
			for jj=1:length(AllGen)
				
				%Find the end node indices
				EndNodeInd=find(ismember(lower(YbusOrderVect),lower(buslist(AllGen(jj)))));												% Find Indices of the individual end node
				
				%Find the phases of end node which match the CN
				EndPhasesInd=find(ismember(YbusPhaseVect(EndNodeInd),YbusPhaseVect(Ind)));
				
				%Find the end nodes indices with PV 
				PvPhases=find(PV(EndNodeInd));
				
				%Match the phases of CN with the corresponding phases of PV
				MatchPVandCN=find(ismember(YbusPhaseVect(Ind),YbusPhaseVect(EndNodeInd(PvPhases))));
				
				
				PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=PvMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+PV(EndNodeInd(EndPhasesInd));			% Update Pv
				LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))=LoadMat(YbusPhaseVect(EndNodeInd(EndPhasesInd)))+LOAD(EndNodeInd(EndPhasesInd));	% update load

				%Update weights
				Weights(Ind(MatchPVandCN),:)=Weights(Ind(MatchPVandCN),:)+Weights(EndNodeInd(PvPhases),:);																		% Update weights
			end
			LOAD(Ind)=LOAD(Ind)+LoadMat(YbusPhaseVect(Ind));
			PV(Ind)=PV(Ind)+PvMat(YbusPhaseVect(Ind));

			%Store rows to delete
			Store=[Store; AllGen];
		end
	end
end

%Update store to reflect Ybus order
Store=unique(Store);
Reduction(2)=length(Store);
S=find(ismember(lower(YbusOrderVect),buslist(Store)));

%Reduce everything to collapse chitlins
Ybus(S,:)=[]; Ybus(:,S)=[];
Lines(find(ismember(lower(Lines(:,2)),buslist(Store))),:)=[];
Lines(find(ismember(lower(Lines(:,3)),buslist(Store))),:)=[];
buslist(Store)=[];
PV(S)=[];
LOAD(S,:)=[];
YbusOrderVect(S)=[];
YbusPhaseVect(S)=[];
Node_number=[];
for ii=1:length(buslist)
Ind=find(strcmpi(buslist(ii),YbusOrderVect))';
Node_number(Ind)=ii;
end
Weights(S,:)=[];


%Need to update critical numbers to busnames
Critical_numbers=find(ismember(buslist,critical_nodes));
Store=[];

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
			Ind1=find(strcmpi(buslist(ParentNodes(ii)),YbusOrderVect))';
			
			%The ind of the critical nodes
			IndCrit=find(strcmpi(buslist(ParentNodes(1)),YbusOrderVect))';
			
			%The ind of the parent node on the other side
			Ind2=find(strcmpi(buslist(ParentNodes(ii+1)),YbusOrderVect))';
			
			%Find out if any missed connections make a matrix of
			%connections work correspondign to that.
			
			%phases        a    b    c
			%         Crit  x    x    x
			%         Ind1  x    x    x
			%         Ind2  x    x    x
			
			ConnMat=zeros(3);
			ConnMat(1,YbusPhaseVect(IndCrit))=IndCrit;  %set row crit nodes
			ConnMat(2,YbusPhaseVect(Ind1))=Ind1;		%set row midd nodes
			ConnMat(3,YbusPhaseVect(Ind2))=Ind2;		%set row next nodes

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
						PV(IndCrit(iii))=PV(IndCrit(iii))+PV(Ind1(iii));
						LOAD(IndCrit(iii),:)=LOAD(IndCrit(iii),:)+LOAD(Ind1(iii),:);
						Weights(IndCrit(iii),:)=Weights(IndCrit(iii),:)+Weights(Ind1(iii),:);
						
					%if middle and next, push everything to next node
					elseif ConnMat(3,1)>0 && ConnMat(2,1)>0
						PV(Ind2(iii))=PV(Ind2(iii))+PV(Ind1(iii));
						LOAD(Ind2(iii))=LOAD(Ind2(iii))+LOAD(Ind1(iii));
						Weights(Ind2(iii),:)=Weights(Ind2(iii),:)+Weights(Ind1(iii),:);
					end
				end
			end
			
			%if full matrix, modify nada
			if sum(ColFlag)==3
				Ind1mod=Ind1;
				Ind2mod=Ind2;
				IndCritmod=IndCrit;
				Z1=inv(Ybus(Ind1,IndCrit));
				Z2=inv(Ybus(Ind1,Ind2));
				
			%else if 1st and 3rd phase connectd
			elseif ~any(ColFlag-[1 0 1])
				Ind1mod=[Ind1(1) Ind1(end)];
				Ind2mod=[Ind2(1) Ind2(end)];
				IndCritmod=[IndCrit(1) IndCrit(end)];
				Z1=inv(Ybus(Ind1mod,IndCritmod));
				Z2=inv(Ybus(Ind1mod,Ind2mod));
				
			%if first and second
			elseif ~any(ColFlag-[1 1 0])
				Ind1mod=[Ind1(1) Ind1(2)];
				Ind2mod=[Ind2(1) Ind2(2)];
				IndCritmod=[IndCrit(1) IndCrit(2)];
				Z1=inv(Ybus(Ind1mod,IndCritmod));
				Z2=inv(Ybus(Ind1mod,Ind2mod));	
				
			% if 2nd and 3rd
			elseif ~any(ColFlag-[0 1 1])
				%ind1
				if length(Ind1)==2
				Ind1mod=[Ind1(1) Ind1(2)];
				else
				Ind1mod=[Ind1(2) Ind1(3)];
				end
				%ind2
				if length(Ind2)==2
				Ind2mod=[Ind2(1) Ind2(2)];
				else
				Ind2mod=[Ind2(2) Ind2(3)];
				end
				%crit
				if length(IndCrit)==2
				IndCritmod=[IndCrit(1) IndCrit(2)];
				else
				IndCritmod=[IndCrit(2) IndCrit(3)];
				end
				Z1=inv(Ybus(Ind1mod,IndCritmod));
				Z2=inv(Ybus(Ind1mod,Ind2mod));
				
			%if only one phase
			else									%If only one phase is connected
				Ind2mod=ConnMat(3,find(ColFlag));
				Ind1mod=ConnMat(2,find(ColFlag));
				IndCritmod=ConnMat(1,find(ColFlag));
				Z1=inv(Ybus(Ind1mod,IndCritmod));
				Z2=inv(Ybus(Ind1mod,Ind2mod));
			end
			
			%get equivalent Z
			Zeq=Z1+Z2;
			
			%rewrite ybus
			Ybus(IndCritmod,Ind2mod)=inv(Zeq);
			Ybus(Ind2mod,IndCritmod)=inv(Zeq);
			
			%calculate weightings
			M2=Z2*inv(Zeq);	M1=Z1*inv(Zeq);

			PV(IndCritmod)=PV(IndCritmod)+(M2*PV(Ind1mod));
			Weights(IndCritmod,:)=Weights(IndCritmod,:)+(M2*Weights(Ind1mod,:));
			LOAD(IndCritmod)=LOAD(IndCritmod)+(M2*LOAD(Ind1mod));
			
			PV(Ind2mod)=PV(Ind2mod)+(M1*PV(Ind1mod));
			Weights(Ind2mod,:)=Weights(Ind2mod,:)+(M1*Weights(Ind1mod,:));
			LOAD(Ind2mod)=LOAD(Ind2mod)+(M1*LOAD(Ind1mod));
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

%Update store to reflect Ybus order
Store=unique(Store);
Reduction(3)=length(Store);
S=find(ismember(lower(YbusOrderVect),buslist(Store)));

%Reduce everything to collapse chitlins
Ybus(S,:)=[]; Ybus(:,S)=[];
buslist(Store)=[];
PV(S)=[];
LOAD(S,:)=[];
YbusOrderVect(S)=[];
YbusPhaseVect(S)=[];
Node_number=[];
for ii=1:length(buslist)
Ind=find(strcmpi(buslist(ii),YbusOrderVect))';
Node_number(Ind)=ii;
end
Weights(S,:)=[];
%Need to update critical numbers to busnames
Critical_numbers=find(ismember(buslist,critical_nodes));
Store=[];

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

tic
fprintf('\n\nFinally rewriting the circuit: ')

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
	circuit.pvsystem(ii).Kv=circuit.circuit.basekv/sqrt(3);
end

%Load
%get PV kVa and pf
BusWithLOAD=find(LOAD>0);
circuit.load=dssload;
% BusListNumWithLOAD=find(ismember(buslist,unique(YbusOrderVect(BusWithLOAD))));
for ii=1:length(BusWithLOAD)
	circuit.load(ii)=dssload;
% 	Ind=find(ismember(YbusOrderVect,buslist(BusListNumWithLOAD(ii))));
% 	PhaseStr=[];
% 	for j=1:length(Ind)
% 		PhaseStr=[PhaseStr '.' num2str(YbusPhaseVect(Ind(j)))];
% 	end
	circuit.load(ii).Name=['LOAD_on_' char(YbusOrderVect(BusWithLOAD(ii))) '.' num2str(YbusPhaseVect(BusWithLOAD(ii)))];%['Load_on_' char(buslist(BusListNumWithLOAD(ii))) PhaseStr];
	circuit.load(ii).phases=1;
	circuit.load(ii).bus1=[char(YbusOrderVect(BusWithLOAD(ii))) '.' num2str(YbusPhaseVect(BusWithLOAD(ii)))];
% 	circuit.load(ii).Pf=num2str(real(LOAD(BusWithLOAD(ii)))./sqrt(real(LOAD(BusWithLOAD(ii))).^2+imag(LOAD(BusWithLOAD(ii))).^2));
	circuit.load(ii).Kw=num2str(real(LOAD(BusWithLOAD(ii))));
	circuit.load(ii).Kvar=num2str(imag(LOAD(BusWithLOAD(ii))));
	circuit.load(ii).Kv=circuit.circuit.basekv/sqrt(3);
end

% %delete useless classes now
names=fieldnames(circuit);
keepFields={'load','buslist','line','circuit','capcontrol','linecode','transformer','capacitor','basevoltages','regcontrol','pvsystem'};
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
		count=count+1;
		Buses{count,1}=busTo;
		Buses{count,2}=BusFrom{jj};
	end
end
Store=zeros(length(Buses),1);
for ii=1:length(Buses)
	if all(ismember([Buses{ii,1} Buses{ii,2}],lower(xfrmrBuses)))
		Store(ii)=1;
	end
end
Buses(find(Store),:)=[];
Store=[];
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

%update line and linecode
circuit.linecode=dsslinecode;
for ii=1:length(Buses)
	
	circuit.linecode(ii)=dsslinecode;
	circuit.linecode(ii).Name=['NewLineCode_bus_' char(Buses{ii,1}) '_to_' char(Buses{ii,2})];
	
	circuit.linecode(ii).R1=[];
	circuit.linecode(ii).R0=[];
	circuit.linecode(ii).X0=[];
	circuit.linecode(ii).X1=[];
	circuit.linecode(ii).Units='kft';
	
	%get r and x matrix
	Ind1= find(ismember(lower(YbusOrderVect),Buses{ii,1}));
	Ind2=find(ismember(lower(YbusOrderVect),Buses{ii,2}));
	
	%Make sure Phases are in correct order
	[~,I]=sort(YbusPhaseVect(Ind1));
	Ind1=Ind1(I);
	[~,I]=sort(YbusPhaseVect(Ind2));
	Ind2=Ind2(I);
	
	Mat=full(Ybus(Ind1,Ind2));
	missingrow=~any(Mat,2);
	missingcol=~any(Mat,1);
	Mat(missingrow,:)=[];
	Mat(:,missingcol)=[];
	FullMat=-inv(Mat);
	
	id1=find(ismember(lower(Lines(:,2)),char(Buses{ii,1})));
	id2=find(ismember(lower(Lines(:,3)),char(Buses{ii,2})));
	FullMat=FullMat./cell2mat(Lines(id1(find(ismember(id1,id2))),1))*1000;
% 	if LineLength(ii)>0
% 	FullMat=FullMat./LineLength(ii)*1000;
% 	end
	
% % % 	%Symmetric definition
% % % 	if length(FullMat)==1
% % % 		circuit.linecode(ii).R1=real(FullMat(1,1));
% % % 		circuit.linecode(ii).X1=imag(FullMat(1,1));
% % % 		circuit.linecode(ii).R0=real(FullMat(1,1))/3;
% % % 		circuit.linecode(ii).X0=imag(FullMat(1,1))/3;
% % % 	else
% % % 		circuit.linecode(ii).R0=real(FullMat(1,1)+2*FullMat(1,2));
% % % 		circuit.linecode(ii).X0=imag(FullMat(1,1)+2*FullMat(1,2));
% % % 		circuit.linecode(ii).R1=real(FullMat(1,1)-FullMat(1,2));
% % % 		circuit.linecode(ii).X1=imag(FullMat(1,1)-FullMat(1,2));
% % % 	end

	%Matrix defintion
	circuit.linecode(ii).NPhases=length(FullMat);
	if length(FullMat)==1
		circuit.linecode(ii).Rmatrix=['(' num2str(real(FullMat(1,1))) ')'];
		circuit.linecode(ii).Xmatrix=['(' num2str(imag(FullMat(1,1))) ')'];
	elseif length(FullMat)==2
		circuit.linecode(ii).Rmatrix=['(' num2str(real(FullMat(1,1))) '|' num2str(real(FullMat(2,1:2)))  ')'];
		circuit.linecode(ii).Xmatrix=['(' num2str(imag(FullMat(1,1))) '|' num2str(imag(FullMat(2,1:2)))  ')'];
	else
		circuit.linecode(ii).Rmatrix=['(' num2str(real(FullMat(1,1))) '|' num2str(real(FullMat(2,1:2))) '|' num2str(real(FullMat(3,1:3))) ')'];
		circuit.linecode(ii).Xmatrix=['(' num2str(imag(FullMat(1,1))) '|' num2str(imag(FullMat(2,1:2))) '|' num2str(imag(FullMat(3,1:3))) ')'];
	end

end

circuit.line=dssline;
for ii=1:length(Buses)
	circuit.line(ii).Name=[char(Buses{ii,1}) '_' char(Buses{ii,2})];
	circuit.line(ii).Linecode=circuit.linecode(ii).Name;
	
	Phases1=YbusPhaseVect(find(ismember(lower(YbusOrderVect),lower(Buses{ii,1}))));
	Phases2=YbusPhaseVect(find(ismember(lower(YbusOrderVect),lower(Buses{ii,2}))));
	Phases=Phases1(find(ismember(Phases1,Phases2)));
	PhaseStr=[];
	for j=1:length(Phases)
		PhaseStr=[PhaseStr '.' num2str(Phases(j))];
	end
	
	circuit.line(ii).bus1=[char(Buses{ii,1}) PhaseStr];
	circuit.line(ii).Phases=length(Phases);
	circuit.line(ii).bus2=[char(Buses{ii,2}) PhaseStr];
	id1=find(ismember(lower(Lines(:,2)),char(Buses{ii,1})));
	id2=find(ismember(lower(Lines(:,3)),char(Buses{ii,2})));
	circuit.line(ii).Length=Lines(id1(find(ismember(id1,id2))),1);
	circuit.line(ii).Units='ft';
	circuit.line(ii).R1=[];
	circuit.line(ii).R0=[];
	circuit.line(ii).X0=[];
	circuit.line(ii).X1=[];
	circuit.line(ii).Rmatrix=[];%circuit.linecode(ii).Rmatrix;
	circuit.line(ii).Xmatrix=[];%circuit.linecode(ii).Xmatrix;
end

%update capcontrol to match updated line names
if isfield(circuit,'capcontrol')
	for ii=1:length(capConBuses1)
		Numms1=find(ismember([Buses{:,1}],char(capConBuses1{ii})));
		Numms2=find(ismember([Buses(:,2)],char(capConBuses2{ii})));
		circuit.capcontrol(ii).Element=['line.' circuit.line(ismember(Numms1,Numms2)).Name];
	end
end

t_=toc;
fprintf('time elapsed %f',t_)

%write weights
circuit.weights=Weights;
circuit.PvbusMap=PVbusMap;
circuit.PvSize=PvSize;
for ii=1:length(YbusOrderVect)
nodeNames{ii}=[YbusOrderVect{ii} '.' num2str(YbusPhaseVect(ii))];
end
circuit.CriticalNode=nodeNames;
circuit.PhaseReduction=Reduction;

outputdss=dsswrite(circuit,[],0,pwd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
end

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
	
	numval = rpncal(val);
	if isempty(numval)
		obj.(prop) = val;
	else
		obj.(prop) = numval;
	end
end
end

function pathtofile = dsswrite( dsscircuit, filename, splitFileFlag, savepath, commands)
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
	savepath = [pwd '/tmp/'];
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