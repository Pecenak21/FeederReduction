function [neighbor] = neighbor(nodeOfInterest, graph)
% nodeOfInterest is a set of nodes (e.g. PV, Storage) that contains all the
% nodes for neighbor finding. All neighbors are in this set.
%
% output: neighborIndex : indices of neighbors based on the nodeOfInterest
% set's index
%%
noi = nodeOfInterest;
g = graph;
%% all nodes' IDs
nid = cell(1,length(g.Nodes));
for i = 1:length(g.Nodes)
	nid{i} = g.Nodes(i).ID;
end
%%
nb.nodeOfInterest = noi;
%% indices of nodes of interest
[v, noiIdx] = ismember(noi,nid);

%% neighbor matrix
nm = zeros(length(noi),length(noi));
weight = zeros(length(noi),length(noi));
for j = 1:length(noi)
	for i = setdiff(1:length(noi),j)
		[v, p] = shortestpath(g,noiIdx(j),noiIdx(i));
		if ~isempty(p)
			% if they are next to each other or if the path doesn't contain
			% any nodes in the set of NOI, then they are neighbors
			[a, b] = ismember(p(2:end-1),noiIdx);
			if length(p) == 2 || (length(p)>2 && ~any(b))
				nm(j,i) = 1;
				weight(j,i) = v;
			end
		end
	end
end
nb.neighborMatrix = nm | eye(size(nm));
nb.weight = weight;
%%
neighbor = nb;
end

