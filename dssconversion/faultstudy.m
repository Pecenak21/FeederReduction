function faultstudy(c1,d1,p1,flgPlot)

if nargin<4
    flgPlot=1;
end

%% setup parameters
cc = c1;
dd = d1;
fs_mode = '3ph'; % 3ph or ll or lg

% The next two sections are experimental
%% Remove regulators?
if(isfield(cc,'regcontrol'))
	% find the transformers
	tid = ismember({cc.transformer.Name},{cc.regcontrol.transformer});
	% replace them with lines
	for i=find(tid);
		nl = dssline('Name',cc.transformer(i).Name);
		[nl.bus1 nl.bus2] = cc.transformer(i).Buses;
		cc.line(end+1) = nl;
	end
	% delete the transformers
	cc.transformer(tid)=[];
	if(isempty(cc.transformer))
		cc = rmfield(cc,'transformer');
	end
	% and delete the regcontrollers
	cc = rmfield(cc,'regcontrol');
end
p = dsswrite(cc,cc.circuit.bus1(1:end),1,['tmp/' cc.circuit.bus1(1:end)]);

%% 
if(isfield(cc,'regcontrol'))
	% find the transformers
	tid = ismember({cc.transformer.Name},{cc.regcontrol.transformer});
	for i=find(tid)
		cc.transformer(i).Conns = {'wye','wye'};
	end
end
cc.generator.Enabled = 'no';

%% load fault data for circuit
switch lower(cc.circuit.Name(1:3))
	case 'cab'
		cnum = '0480';
	case 'alp'
		cnum = '0355';
	case 'avo'
		cnum = '0520';
	case 'val'
		cnum = '0909';
	case 'cre'
		cnum = '0971';
	otherwise 
		error;
end
x = excel2obj(sprintf('dssconversion/custdata/27_%s_-_Sections.csv', cnum) );
x = struct2cell(x); x = x{1};

% convert section names to the format we're using for the DSS objects
x = structconv(x);
x.Section_Id = regexprep(x.Section_Id,' ','_');
fns = {'Symmetrical_Amps_LG_Min','Symmetrical_Amps_LG_Max','Symmetrical_Amps_LL','Symmetrical_Amps_LLG','Symmetrical_Amps_3Ph','Asymmetrical_Amps_LL'};
for i = fns; i = i{1};
	a = cellfun(@ischar,x.(i));
	[x.(i){a}] = deal(NaN);
end
x = structconv(x);

%% load the fault data calculated in opendss
curDir = pwd;
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
pause(1);
z = faultread(['tmp/' subdir '/' cc.circuit.Name '_FaultStudy.Txt']);
cd(curDir);
%% Map Utility data to busid
[a b] = ismember({x.Section_Id},{cc.line.Name});
if(~all(a))
	warning('faults:missing',['Some sectionIds are present in the fault data that are missing in the network data:' sprintf('\n\t%s',x(~a).Section_Id)]);
	b = b(a);
end
l = regexprep(cc.line(b).bus2,'(\.\d+)+$','');
[x(a).BusID] = deal(l{:});
l = regexp({x(~a).Section_Id}','[^_]+$','match','once');
[x(~a).BusID] = deal(l{:});

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
[b c] = ismember({x.BusID},strtrim({y.bus}'));
if(~all(b))
	warning('faults:missing',['Some busIds are present in lines for the fault data, but not in OpenDSS''s fault output' sprintf('\n\t%s',x(~b).BusID)]);
	c = c(b);
end
% extract the current data; column 1 will be synergee, column 2 opendss
switch(lower(fs_mode))
	case 'll'
		current = [vertcat(x(b).Symmetrical_Amps_LL) vertcat(y(c).LL)];
	case '3ph'
		current = [vertcat(x(b).Symmetrical_Amps_3Ph) vertcat(y(c).I3Phase)];
	case {'lg','1ph'}
		current = [vertcat(x(b).Symmetrical_Amps_LG_Max) vertcat(y(c).I1Phase)];
end
dist = vertcat(x(b).Dist_kFt);
% lookup xy coordinates
[xy xy] = ismember(strtrim({y(c).bus}),cc.buslist.id);
xy = cc.buslist.coord(xy,:);
% and calculate the fractional difference
cdiff = current(:,2)./current(:,1)-1;
% and get a sort order (now named b) for sorting by distance
[b b] = sort(dist);
%% Plot
figure;
plot(dist(b),current(b,1),'x',dist(b),current(b,2),'.');
legend('Source Data, Symmetrical Amps LL','OpenDSS');
xlabel('Distance, kft')
ylabel('Current, A');
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
legend(h(2:3),{'Source Data, Symmetrical Amps LL','OpenDSS'});
xlabel('Distance, kft')
ylabel('Current, A');
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
[ism_ locm] = ismember(m,strtrim({y(c).bus}));
[isn_ locn] = ismember(n,strtrim({y(c).bus}));
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

%% Manual fault calculations at the substation:
% clear i;
me(1) = dd.InstFeeders.NominalKvll/sqrt(3)/(dd.InstFeeders.PosSequenceResistance+j*dd.InstFeeders.PosSequenceReactance);
me(2) = dd.InstFeeders.NominalKvll/sqrt(3)*3/(dd.InstFeeders.PosSequenceResistance*2+dd.InstFeeders.ZeroSequenceResistance+j*(dd.InstFeeders.PosSequenceReactance*2+dd.InstFeeders.ZeroSequenceReactance));
me(3) = dd.InstFeeders.NominalKvll/(dd.InstFeeders.PosSequenceResistance*2+j*(dd.InstFeeders.PosSequenceReactance*2));
abs(me);
