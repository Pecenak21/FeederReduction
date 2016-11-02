function faultstudy_Vadim(c1,d1,p1,flgPlot)

if nargin<4
    flgPlot=1;
end

%% setup parameters
cc = c1;
dd = d1;
p = p1;
fs_mode = '3ph'; % 3ph or ll or lg

% The next two sections are experimental
%% Remove regulators?
% if(isfield(cc,'regcontrol'))
% 	% find the transformers
% 	tid = ismember({cc.transformer.Name},{cc.regcontrol.transformer});
% 	% replace them with lines
% 	for i=find(tid);
% 		nl = dssline('Name',cc.transformer(i).Name);
% 		[nl.bus1 nl.bus2] = cc.transformer(i).Buses;
% 		cc.line(end+1) = nl;
% 	end
% 	% delete the transformers
% 	cc.transformer(tid)=[];
% 	if(isempty(cc.transformer))
% 		cc = rmfield(cc,'transformer');
% 	end
% 	% and delete the regcontrollers
% 	cc = rmfield(cc,'regcontrol');
% end
% p = dsswrite(cc,cc.circuit.bus1(2:end),1,cc.circuit.bus1(2:end));

%% 
% % if(isfield(cc,'regcontrol'))
% % % 	find the transformers
% % 	tid = ismember({cc.transformer.Name},{cc.regcontrol.transformer});
% % 	for i=find(tid)
% % 		cc.transformer(i).Conns = {'wye','wye'};
% % 	end
% % end
% % cc.generator.Enabled = 'no';

%% load fault data for circuit
[subdir subdir] = fileparts(fileparts(p));
if strcmp(subdir,'SCE_Centaur')
    % x = excel2obj(sprintf('custdata/27_%s_-_Sections.csv', cc.circuit.bus1));
    x = excel2obj(sprintf('SC-Centaur-SPVP022-2013-01-30.xlsx'));
elseif strcmp(subdir,'SCE_Durox')
    x = excel2obj(sprintf('Durox_SCDuty_130613_Gen.xlsx'));  
else
    warning('SC File not found');
    return
end

x = struct2cell(x); x = x{1};

% convert section names to the format we're using for the DSS objects
x = structconv(x);
% x.Section_Id = regexprep(x.Section_Id,' ','_');
x.Node_Id_ = regexprep(x.Node_Id_,' ','_');
% fns = {'Symmetrical_Amps_LG_Min','Symmetrical_Amps_LG_Max','Symmetrical_Amps_LL','Symmetrical_Amps_LLG','Symmetrical_Amps_3Ph','Asymmetrical_Amps_LL'};
fns = {'LG_Amps_','LL_Amps_','LLG_Amps_','LLL_Amps_'};
for i = fns; i = i{1};
	a = cellfun(@ischar,x.(i));
	[x.(i){a}] = deal(NaN);
end
x = structconv(x);

%% load the fault data calculated in opendss
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
% UtilityBus1 = regexprep({cc.line.bus1},'(\.\d+)+$','');
% UtilityBus2 = regexprep({cc.line.bus2},'(\.\d+)+$','');
% [a b] = ismember({x.Node_Id_},UtilityBus1);
% if(~all(a))
% 	warning('faults:missing',['Some sectionIds are present in the fault data that are missing in the network data:' sprintf('\n\t%s',x(~a).Node_Id_)]);
% 	b = b(a);
% end
% % l = regexprep(UtilityBus1(b),'(\.\d+)+$','');
% % [x(a).BusID] = deal(l{:});
% [x(a).BusID] = deal(UtilityBus1{b});
% l = regexp({x(~a).Node_Id_}','[^_]+$','match','once');
% [x(~a).BusID] = deal(l{:});

%% match buses
% calculate a mask (b) for the synergee data, and a lookup for matching
% busses (c(b)) in the opendss data
% [b c] = ismember({x.BusID},{y.bus}');
% x.Node_Id_ = regexprep(x.Node_Id_,' ','_');
% x.Node_Id_ = regexprep(x.Node_Id_,'\$','_');
NOD = regexprep({x.Node_Id_},'\-','_');
[b c] = ismember(NOD,{y.bus}');
if(~all(b))
	warning('faults:missing',['Some busIds are present in lines for the fault data, but not in OpenDSS''s fault output' sprintf('\n\t%s',x(~b).Node_Id_)]);
	c = c(b);
end
% extract the current data; column 1 will be synergee, column 2 opendss
switch(lower(fs_mode))
	case 'll'
		current = [vertcat(x(b).LL_Amps_) vertcat(y(c).LL)];
	case '3ph'
% 		current = [vertcat(x(b).LLL_Amps_) vertcat(y(c).I3Phase)];
        current = [vertcat(x(b).LLG_Amps_) vertcat(y(c).I3Phase)];
	case {'lg','1ph'}
		current = [vertcat(x(b).LG_Amps_) vertcat(y(c).I1Phase)];
end
LinesUtility = {x(b).Node_Id_}';
LinesOpenDSS = {y(c).bus}';

dist = vertcat(x(b).Total_distance_ft);
% lookup xy coordinates
[xy xy] = ismember({y(c).bus},cc.buslist.id);
% xy = cc.buslist.coord(xy,:);
% and calculate the fractional difference
% cdiff = current(:,2)./current(:,1)-1;
% and get a sort order (now named b) for sorting by distance
[b b] = sort(dist);
%% Plot
% figure;
% plot(dist(b),current(b,1),'x',dist(b),current(b,2),'.');
% legend('Source Data, Symmetrical Amps LG','OpenDSS');
% xlabel('Distance, kft')
% ylabel('Current, A');
%% Plot Current vs distance from substation
% Arrange the data so that we get disconnected vertical lines:
% we do this by using the (x,y1,x,y2,nan,nan) for each point, and then
% reshaping the matrix.
dat_ = [dist(b) current(b,1) dist(b) current(b,2) nan(length(dist),2)];
dat_ = reshape(dat_',2,numel(dat_)/2)';
% plot the data in a new figure
figure;
h = plot(dat_(:,1)/1000, dat_(:,2),'r',dist(b)/1000,current(b,1),'x',dist(b)/1000,current(b,2),'.');
TEST = {};
test = {LinesUtility(b),LinesOpenDSS(b),current(b,1),current(b,2)};
TEST = [test;TEST];

% set some labels and make sure we don't overwrite this plot
legend(h(2:3),{'Source Data','OpenDSS'});
xlabel('Distance, kft')
ylabel('Current, A');
set(gcf,'nextplot','new')
grid on;
%%
%save fault1.mat dist current xy
%%
% plot(dist(b),cdiff(b),'x');
%%
% mask = ~isnan(cdiff);
% figure;
% plot3(xy(mask,1),xy(mask,2),cdiff(mask),'x');

%% Draw the circuit in 3D with difference on the z-axis
% lookup lines and map the data onto them
% m = regexprep({cc.line.bus1}','(\.\d+)+$','');
% n = regexprep({cc.line.bus2}','(\.\d+)+$','');
% [ism_ locm] = ismember(m,{y(c).bus});
% [isn_ locn] = ismember(n,{y(c).bus});
% mask = ism_ & isn_;
% locm = locm(mask);
% locn = locn(mask);
% collect the data for plotting
% (x1,y1,z1,x2,y2,z2,nan,nan,nan) for each line, then reshape to be able to
% plot as a single object
% dat_ = [xy(locm,:),cdiff(locm),xy(locn,:),cdiff(locn),nan(length(locm),3)];
% dat_ = reshape(dat_',3,numel(dat_)/3)';
% figure;
% plot3(dat_(:,1),dat_(:,2),dat_(:,3))
% xlabel('X'); ylabel('Y'); zlabel('pu \Delta I_{sc}');
% set(gca,'looseinset',get(gca,'tightinset'));

%% Manual fault calculations at the substation:
% clear i;
% me(1) = dd.InstFeeders.NominalKvll/sqrt(3)/(dd.InstFeeders.PosSequenceResistance+j*dd.InstFeeders.PosSequenceReactance);
% me(2) = dd.InstFeeders.NominalKvll/sqrt(3)*3/(dd.InstFeeders.PosSequenceResistance*2+dd.InstFeeders.ZeroSequenceResistance+j*(dd.InstFeeders.PosSequenceReactance*2+dd.InstFeeders.ZeroSequenceReactance));
% me(3) = dd.InstFeeders.NominalKvll/(dd.InstFeeders.PosSequenceResistance*2+j*(dd.InstFeeders.PosSequenceReactance*2));
% abs(me);
