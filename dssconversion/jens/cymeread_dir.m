function [ ob ] = cymeread_dir( prefix )
% Convert cyme files with multiple sheets to object oriented struct
fs = dir([prefix '*']);
if(length(fs)==1 && isdir(prefix))
	fs = dir([prefix '/']);
else
	prefix = fileparts(prefix);
end
fns = strcat(prefix,'/',{fs.name}');
fns(cellfun(@isempty,regexpi(fns,'(?:equipment|load|network).txt$','once'))) = [];

ob = cymeread(fns);

% do some post-processing/cleanup
ob = sanitizeIds(ob);
tic;
ob = do_net_ids(ob);
toc
return;
ob = do_net_ids(sanitizeIds(ob));

end

function ob = sanitizeIds(ob)
	fns = fieldnames(ob);
	for i=1:length(fns)
        if(length(ob) && isstruct(ob(1).(fns{i})))
			ob.(fns{i}) = sanitizeIds(ob.(fns{i}));
        elseif(length(ob) && iscell(ob(1).(fns{i})))
            for i_=1:length(ob(1).(fns{i}))
                ob.(fns{i}){i_} = sanitizeIds(ob.(fns{i}){i_});
            end
		elseif(regexpi(fns{i},'(section|node)id'))
			n = length(unique({ob.(fns{i})}));
			x = fnSanitize({ob.(fns{i})});
			[ob.(fns{i})] = deal(x{:});
			if(length(unique(x)) ~= n)
				warning('ID sanitization may have caused id conflicts!');
			end
		end
	end
end

function ob = do_net_ids(ob)
	try
		heads = ob.network.headnodes;
	catch
		try
			heads = ob.network.section.feeder;
		catch
			return;
		end
		heads.NodeID = deal(heads.HeadNodeID);
	end

	f = {ob.network.section.section.FromNodeID}';
	t = {ob.network.section.section.ToNodeID}';
	[tfmap, tfmap] = ismember(f,t); % tfmap now becomes a map of which section is connected upstream of a given section.
	% our strategy is now to start for each node and work up the chain until we find either a node that has a NetworkID assigned, or one that has no upstream nodes (in which case we hope that node is a head node and we look up a NetworkID for it in the appropriate table
	nid = cell(size(f));
    nid(:)={''}; % js: added this line to ensure that Network IDs are strings only, mixed string/numeric arrays create problems later on
	badheads = {};
	% start by assigning the correct network id for each head node
	[i,i] = ismember({heads.NodeID},f);
	if(~all(i)), warning('not all head node ids exist in the network'); end
	nid(i(i~=0)) = {heads(i~=0).NetworkID};

	for i=1:length(f)
		if(~isempty(nid{i})); continue; end
		j = i;
		while(tfmap(j(1)) && isempty(nid{j(1)}) && ~any(j(1)==j(2:end)))
			j = [tfmap(j(1)) j];
		end
		% for the case where we hit a top, try to look up this as a head node
		if(~tfmap(j(1)))
			k = find(strcmp(f{j(1)},{heads.NodeID}));
			if(length(k)==0)
				if(~any(strcmp(f{j(1)},badheads)))
					warning('found a chain head with no network ID! "%s"',f{j(1)})
					badheads = [badheads f{j(1)}];
				end
				continue;
			end
			k = heads(k(1)).NetworkID;
		elseif(~isempty(nid{j(1)})) % found one with a network id already
			k = nid{j(1)};
			j = j(2:end);
		else
			warning('may have found infinite loop; unclear how this will work');
			% try to remove the loop
			[xx,xx] = unique(t);
			tt = t; tt(xx) = [];
			tt = find(ismember(t,tt));
			if(length(tt)==2 && sum(ismember(tt,j))==1) % some conditions in which this tweak works well; otherwise we'll have other issues probably
				yy = find(tfmap==tt(ismember(tt,j)));
				if(sum(ismember(yy,j))==1)
					tfmap(yy(ismember(yy,j))) = tt(~ismember(tt,j));
					warning('tried to remove the loop');
				end
			end
			continue;
		end
		% set these nodes to have the found network ID
		nid(j) = {k};
	end

	fprintf('NetworkID found for %i of %i nodes',sum(~cellfun(@isempty,nid)),length(nid));
	% assign the data into our structure
	[ob.network.section.section.NetworkID] = deal(nid{:});
end

