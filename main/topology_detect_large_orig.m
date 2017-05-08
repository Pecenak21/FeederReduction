function [topo, generation]=topology_detect_large_orig(topo,generation,Ybus,parent,Node_number)
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
% column 6: Physical distance to the substation, should be in kilometers

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
adj_bus(find(adj_bus==topo(parent,2)))=0;
% for i=1:max(Node_number)
% 	adj_bus(find(adj_bus==topo(i,2)))=0;
% end
adj_bus(find(adj_bus==0))=[]; %remove parent

% if any(ismember(adj_bus,generation{parent,4}))
% 	error(['Mesh detected on bus ' num2str(parent)]);
% end
generation{parent,2}=[];	  %set children empty
generation{parent,3}=[];	  %set grandchildren empty
if length(adj_bus~=0)
	for k=1:length(adj_bus)
		child=adj_bus(k);
		if child==447
			stop=1;
		end
		if max(max(topo(:,1:2)==child))==0
			topo(child,1)=child;
			topo(child,2)=parent;
			generation{child,1}=child;
			generation{child,4}=[parent;generation{parent,4}];
			generation{child,5}= generation{parent,5}+1;
			[topo, generation]=topology_detect_large_orig(topo,generation,Ybus,child,Node_number);
			topo(parent,3)=topo(parent,3)+1;
			topo(parent,4)=topo(parent,4)+topo(child,4)+1;
			generation{parent,2}=[generation{parent,2};child];
			generation{parent,3}=[generation{parent,3};child;generation{child,3}];
		end
	end
end
end