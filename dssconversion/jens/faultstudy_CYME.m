function faultstudy_CYME(c1,d1,p1,flgPlot)

if nargin<4
    flgPlot=1;
end

%% setup parameters
cc = c1;
dd = d1;
p = p1;
fs_mode = '3ph'; % 3ph or ll or lg

%% load fault data for circuit
x = excel2obj(['custdata\' 'HOL Short Circuit Results_sorted.xlsx']);
x = struct2cell(x); x = x{1};

% convert section names to the format we're using for the DSS objects
x = structconv(x);
x= StandardizeFieldNames(x,1);
x.NodeID= UnifyFormat(x,'NodeID','str');
x.NodeID=fnSanitize(x.NodeID);

% x.Section_Id = Node(x.Section_Id,' ','_');
% fns = {'Symmetrical_Amps_LG_Min','Symmetrical_Amps_LG_Max','Symmetrical_Amps_LL','Symmetrical_Amps_LLG','Symmetrical_Amps_3Ph','Asymmetrical_Amps_LL'};

fns = {'Symmetrical_Amps_LG_Min','Symmetrical_Amps_LG_Max','Symmetrical_Amps_LL','Symmetrical_Amps_LLG','Symmetrical_Amps_3Ph'};
for i = fns; i = i{1};
	a = cellfun(@ischar,x.(i));
	[x.(i){a}] = deal(NaN);
end
x = structconv(x);

%% load the fault data calculated in opendss
[subdir subdir] = fileparts(fileparts(p));
t = dssget(p);
t = t.Text;
t.Command = 'get editor'; olded = t.Result; if(isempty(olded)), olded = 'notepad.exe'; end
t.Command = 'set editor="where.exe"'; % silence output
t.Command = ['cd "' fileparts(p) '"/'];
t.Command = 'Set mode=faultstudy';
t.Command = 'Solve';
t.Command = ['Export faultstudy ' lower(cc.circuit.Name) '_fault.csv'];
y = faultread(t.Result);
t.Command = 'show faults';
t.command = ['cd "' pwd '"'];
t.Command = ['set editor="' olded '"']; % unsilence output
z = faultread([subdir '/' cc.circuit.Name '_FaultStudy.Txt']);

%% Map Utility data to busid
% [a b] = ismember({x.Section_Id},{cc.line.Name});
% if(~all(a))
% 	warning('faults:missing',['Some sectionIds are present in the fault data that are missing in the network data:' sprintf('\n\t%s',x(~a).Section_Id)]);
% 	b = b(a);
% end
% l = regexprep(cc.line(b).bus2,'(\.\d+)+$','');
% [x(a).BusID] = deal(l{:});
% l = regexp({x(~a).Section_Id}','[^_]+$','match','once');
% [x(~a).BusID] = deal(l{:});

[a b] = ismember(cc.buslist.id,{x.NodeID});

if(~all(a))
	warning('faults:missing',['Some Node IDs are present in the fault data that are missing in the network data:' sprintf('\n\t%s',x(~a).NodeID)]);
	b = b(a);
end
% l = regexprep(cc.line(b).bus2,'(\.\d+)+$','');
% [x(a).BusID] = deal(l{:});
% l = regexp({x(~a).NodeID}','[^_]+$','match','once');
% [x(~a).BusID] = deal(l{:});
% 
% % transformers are handled differently: OpenDSS calculates at nodes, and so
% % if we continue to use the downstream node for the Synergee data, the
% % results will be off by quite a bit.  Instead, for sectionIDs in the
% % synergee data that would match a transformer, we replace that data point
% % with the current for the line connected to the secondary of the
% % transformer.  (We hope that there is only one such line!)
% for i=1:length(cc.transformer)
% 	bid = regexprep(cc.transformer(i).Buses{2},'(\.\d+)+$','');
% 	bidx = find(strcmp(bid,{x.BusID}'));
% 	lidx = find(strcmp(bid,regexprep({cc.line.bus1}','(\.\d+)+$','')));
% 	if(length(lidx)==1 && length(bidx)==1)
% 		lidx = find(strcmp(regexprep(cc.line(lidx).bus2,'(\.\d+)+$',''),{x.BusID}'));
% 		x(bidx).Symmetrical_Amps_LG_Max = x(lidx).Symmetrical_Amps_LG_Max;
% 		x(bidx).Symmetrical_Amps_LL = x(lidx).Symmetrical_Amps_LL;
% 		x(bidx).Symmetrical_Amps_3Ph = x(lidx).Symmetrical_Amps_3Ph;
% 	else
% 		warning('faultstudy:transformerRemap','remap values for synergee data point %i failed!',bidx);
% 	end
% end

%% match buses
% calculate a mask (b) for the synergee data, and a lookup for matching
% busses (c(b)) in the opendss data
[b c] = ismember({y.bus}',{x.NodeID});
if(~all(b))
	warning('faults:missing',['Some busIds are present in lines for the fault data, but not in OpenDSS''s fault output' sprintf('\n\t%s',x(~b).NodeID)]);
	c = c(b);
end
% extract the current data; column 1 will be synergee, column 2 opendss
switch(lower(fs_mode))
	case 'll'
		current = [vertcat(x(c).Symmetrical_Amps_LL) vertcat(y(b).LL)];
        strTitle='Line-to-Line Fault';
	case '3ph'
		current = [vertcat(x(c).Symmetrical_Amps_3Ph) vertcat(y(b).I3Phase)];
        strTitle='Three-Phase Fault';
	case {'lg','1ph'}
		current = [vertcat(x(c).Symmetrical_Amps_LG_Max) vertcat(y(b).I1Phase)];
        strTitle='Line-to-Ground Fault';
end
dist = vertcat(x(c).Distance_ft);
% lookup xy coordinates
[xy xy] = ismember({y(b).bus},cc.buslist.id);
xy = cc.buslist.coord(xy,:);
% and calculate the fractional difference
cdiff = current(:,2)./current(:,1)-1;
% and get a sort order (now named b) for sorting by distance
[b b] = sort(dist);
%% Plot

if flgPlot
    figure;
    plot(dist(b),current(b,1),'x',dist(b),current(b,2),'.');
    legend('Source Data, Symmetrical Amps LL','OpenDSS');
    xlabel('Distance, kft')
    ylabel('Current, A');
    title('strTitle');
    %% Plot Current vs distance from substation
    % Arrange the data so that we get disconnected vertical lines:
    % we do this by using the (x,y1,x,y2,nan,nan) for each point, and then
    % reshaping the matrix.
    dat_ = [dist(b) current(b,1) dist(b) current(b,2) nan(length(dist),2)];
    dat_ = reshape(dat_',2,numel(dat_)/2)';
    % plot the data in a new figure
    figure;
    h = plot(dat_(:,1), dat_(:,2),'r',dist(b),current(b,1),'x',dist(b),current(b,2),'.');
    % set some labels and make sure we don't overwrite this plot
    legend(h(2:3),{'Source Data','OpenDSS'});
    xlabel('Distance, kft')
    ylabel('Current, A');
    title('strTitle');
    set(gcf,'nextplot','new')
    %%
    %save fault1.mat dist current xy
    %%
    plot(dist(b),cdiff(b),'x');
    %%
    mask = ~isnan(cdiff);
    figure;
    plot3(xy(mask,1),xy(mask,2),cdiff(mask),'x');

    %% Draw the circuit in 3D with difference on the z-axis
    % lookup lines and map the data onto them
    m = regexprep({cc.line.bus1}','(\.\d+)+$','');
    n = regexprep({cc.line.bus2}','(\.\d+)+$','');
    [ism_ locm] = ismember(m,{y(b).bus});
    [isn_ locn] = ismember(n,{y(b).bus});
    mask = ism_ & isn_;
    locm = locm(mask);
    locn = locn(mask);
    % collect the data for plotting
    % (x1,y1,z1,x2,y2,z2,nan,nan,nan) for each line, then reshape to be able to
    % plot as a single object
    dat_ = [xy(locm,:),cdiff(locm),xy(locn,:),cdiff(locn),nan(length(locm),3)];
    dat_ = reshape(dat_',3,numel(dat_)/3)';
    figure;
    plot3(dat_(:,1),dat_(:,2),dat_(:,3))
    xlabel('X'); ylabel('Y'); zlabel('pu \Delta I_{sc}');
    set(gca,'looseinset',get(gca,'tightinset'));
end

%% Manual fault calculations at the substation:
% clear i;
% me(1) = dd.InstFeeders.NominalKvll/sqrt(3)/(dd.InstFeeders.PosSequenceResistance+i*dd.InstFeeders.PosSequenceReactance);
% me(2) = dd.InstFeeders.NominalKvll/sqrt(3)*3/(dd.InstFeeders.PosSequenceResistance*2+dd.InstFeeders.ZeroSequenceResistance+i*(dd.InstFeeders.PosSequenceReactance*2+dd.InstFeeders.ZeroSequenceReactance));
% me(3) = dd.InstFeeders.NominalKvll/(dd.InstFeeders.PosSequenceResistance*2+i*(dd.InstFeeders.PosSequenceReactance*2));
me(1) = c1.circuit.basekv/sqrt(3)/(c1.circuit.R1+j*c1.circuit.X1);
me(2) = c1.circuit.basekv/sqrt(3)*3/(c1.circuit.R1*2+c1.circuit.R0+j*(c1.circuit.X1*2+c1.circuit.X0));
me(3) = c1.circuit.basekv/(c1.circuit.R1*2+j*(c1.circuit.X1*2));

abs(me)
