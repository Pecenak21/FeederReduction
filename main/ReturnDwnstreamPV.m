function [pvNamesDown, pvNamesUp, LoadNamesDown, LoadNamesUp] = ReturnDwnstreamPV(Bus_upstream,circuit)

%Created by Zachary K. Pecenak on 6/18/2016

%Example input
% load('c:\users\zactus\gridIntegration\results\ValleyCenter_wpv_existing.mat')
% Bus_upstream={'03551325'};
% or Bus_upstream=c.buslist.id(ceil(length(c.buslist.id)*rand(1,1))) this is a random bus
% [c] = FeederReduction(Bus_upstream,c);

%Check both inputs are met
if nargin<2
	error('This requires two inputs, the desired nodes to keep and the feeder circuit')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% get Ybus from OpenDSS for feeder.
fprintf('\nGetting YBUS and buslist from OpenDSS: ')
tic
%remove load and PV from Ybus
circuit_woPV=rmfield(circuit,'pvsystem');
circuit_woPV=rmfield(circuit_woPV,'load');
%load the circuit and generate the YBUS
p = dsswrite(circuit_woPV,[],0,[]); o = actxserver('OpendssEngine.dss');
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

Bus_upstream(find(~ismember(lower(Bus_upstream),lower(buslist))))=[];
if ~any(ismember(lower(Bus_upstream),lower(buslist)))
 	error('The selected bus in the circuit: \n%s', Bus_upstream{~ismember(Bus_upstream,buslist)})
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
	
	%Get the bus and nodes of the PV
	PVbusInd(j)=find(ismember(lower(buslist),lower(PVname{1})));
	
		%Map PV to appropriate buses for later conversion
	PVbusMap(j,PVbusInd(j))=1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get PV data you want
idx=find(ismember(buslist,lower(Bus_upstream)));
[PVDown,c]=find(PVbusMap(:,generation{idx,3}));
if isempty(PVDown)
	fprintf('\n\nThere is no PV Downstream silly fufu :)\n\n\n')
	pvNamesDown=[];
else
pvNamesDown=circuit.pvsystem(PVDown).name;
end

[PVUp,c]=find(PVbusMap(:,generation{idx,4}));
if isempty(PVUp)
	fprintf('\n\nThere is no PV upstream silly fufu :)\n\n\n')
	pvNamesUp=[];
else
pvNamesUp=circuit.pvsystem(PVUp).name;
end

%% Get Load you want
K=regexp([circuit.load(:).bus1],'\.','split','once');
k=[K{:}]; k(find(cellfun('length',k)==1))=[]; 
[LoadDown]=find(ismember(lower(k),lower(buslist(generation{idx,3}))));
if isempty(LoadDown)
	fprintf('\n\nThere is no Load Downstream silly fufu :)\n\n\n')
	LoadNamesDown=[];
else
LoadNamesDown=circuit.load(LoadDown).name;
end

[LoadUp]=find(ismember(lower(k),lower(buslist(generation{idx,4}))));
if isempty(LoadUp)
	fprintf('\n\nThere is no Load Upstream silly fufu :)\n\n\n')
	LoadNamesUp=[];
else
LoadNamesUp=circuit.load(LoadUp).name;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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