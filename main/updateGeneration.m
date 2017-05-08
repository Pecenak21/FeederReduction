function [generation] = updateGeneration(generation,buslist, buslist_org)
%New TOPO detect after reduction

RemoveFromTopo=find(~ismember(buslist_org,buslist));
generation(RemoveFromTopo,:)=[];

keySet=find(ismember(buslist_org,buslist));
%map buslist
% keySet = lower(buslist);
mapObj = containers.Map(keySet,1:length(keySet));


for ii=1:length(generation)
	generation{ii,1}=ii;
	remainingInd=find(~ismember(cell2mat(generation(ii,2)),RemoveFromTopo));
	generation{ii,2}=generation{ii,2}(remainingInd);
	if ~isempty(generation{ii,2})
		generation{ii,2}=cell2mat(values(mapObj,num2cell(generation{ii,2})));
	end
	
	remainingInd=find(~ismember(cell2mat(generation(ii,3)),RemoveFromTopo));
	generation{ii,3}=generation{ii,3}(remainingInd);
	if ~isempty(generation{ii,3})
		generation{ii,3}=cell2mat(values(mapObj,num2cell(generation{ii,3})));
	end
	
	remainingInd=find(~ismember(cell2mat(generation(ii,4)),RemoveFromTopo));
	generation{ii,4}=generation{ii,4}(remainingInd);
	if ~isempty(generation{ii,4})
		generation{ii,4}=cell2mat(values(mapObj,num2cell(generation{ii,4})));
	end
	generation{ii,5}=length(remainingInd);
end

emptyNot=find(~cellfun(@isempty,generation(:,3)));

for ii=1:length(emptyNot)
	Lengths=[generation{generation{emptyNot(ii),3},5}];
	generation{emptyNot(ii),2}=generation{emptyNot(ii),3}(find(Lengths==min(Lengths)));
end
end