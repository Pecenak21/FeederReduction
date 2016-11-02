function ymat = calculateYMatrix(nodesOfInterest, circuit)
% nodesOfInterest: buses of interest (cell of string)
noi = nodesOfInterest;
nn = length(noi); % number of nodes
ct = circuit;

%% set all loads to 0 power
% for i = 1:length(c.load)
% 	c.load(i).Kw = 0;
% 	c.load(i).Kvar = 0;
% end

%% or just remove them all from the circuit
% c = rmfield(c,'load');

%% set source voltage to 1V
vol = .001;
ct.basevoltages = vol;
% change source
% ct.circuit.Phases = 1;
% ct.circuit.bus1 = 'sourcebus.1';
ct.circuit.basekv = vol;

%% stick fault objects to all nodes of interest
f(nn) = dssfault;
for j = 1:nn
    f(j).Name = ['f_' sprintf('%d',j)];
    f(j).bus1 = noi{j};
end
ct.fault = f;

%% initialize 'gmat'
ymat = zeros(nn,nn);

%% Remove regulators?
if(isfield(ct,'regcontrol'))
	% find the transformers
	tid0 = ismember(lower({ct.transformer.Name}),lower({ct.regcontrol.transformer}));
	% replace them with lines
	for i=find(tid0);
		nl = dssline('Name',ct.transformer(i).Name);
		[nl.bus1, nl.bus2] = ct.transformer(i).Buses;
		nl.length = 1;
		nl.units = 'ft';
		ct.line(end+1) = nl;
	end
	% delete the transformers
	ct.transformer(tid0)=[];
	if(isempty(ct.transformer))
		ct = rmfield(ct,'transformer');
	end
	% and delete the regcontrollers
	ct = rmfield(ct,'regcontrol');
    % remove energymeter that is connected to the transformer
    if isfield(ct,'energymeter')
        for i = 1:length(ct.energymeter)
            if strfind(lower(ct.energymeter(i).element),'transformer')
                ct.energymeter(i) = [];
            end
        end
    end
end

%% find transformer that connects to sourcebus
if isfield(ct,'transformer')
	[~,tid] = ismember('sourcebus',lower(cleanBus([ct.transformer.Bus])));
	if tid > 0, 
		tid = floor(tid/2); 
		[~,regid] = ismember(ct.transformer(tid).Name,{ct.regcontrol.transformer});
	end
	
	for k = 1:length(ct.transformer)
		ct.transformer(k).kVs = [vol vol];
    end
end

%% 
if(isfield(ct,'regcontrol'))
	% find the transformers
	tid = ismember({ct.transformer.Name},{ct.regcontrol.transformer});
	for i=find(tid)
		ct.transformer(i).Conns = {'wye','wye'};
	end
end

%% initialize dssengine
o = dsscompile(ct);

%% process at each node
for i = 1:nn
	%% 
    if exist('tid','var') && tid > 0
        % connect node of interest to sub tranx
        setval(o,ct.transformer(tid),'Buses',['[sourcebus ' noi{i}]);
		% reset regulator bus
		if regid>0 && ~isempty(ct.regcontrol(regid).bus)
            setval(o,ct.regcontrol(regid),'bus',noi{i});
		end
    else
        setval(o,'vsource.source','bus1',noi{i});
    end
    
    [~, id] = ismember(noi{i}, {ct.fault.bus1});
    %% get mutual admittance 
    setval(o,ct.fault(id),'enabled','false');
    o.ActiveCircuit.Solution.Solve;
    v = dssgetval(o,ct.fault,'currents');
	%v(abs(v)<.001) = 0;
    ymat(:,i) = complex(v(:,1),v(:,2));
    
    setval(o,ct.fault(id),'enabled','true');
    %% get self admittance
%     o.ActiveCircuit.Solution.Solve;
%     v = dssgetval(o,ct.fault,'currents');
% 	%v(abs(v)<.001) = 0;
%     ymat(i,i) = complex(v(i,1),v(i,2));
end

%% last cleanup
% for i = 1:nn
%     id = find(ymat(i,:)>0);
%     if length(id)>1
%         ymat(i,id) = 0;
%         ymat(i,i) = abs(sum(ymat(i,:)));
%     end
% end

% clean up: Assume current smaller than 10^-6 is 0
%gmat(abs(gmat)<10^-6) = 0;

end