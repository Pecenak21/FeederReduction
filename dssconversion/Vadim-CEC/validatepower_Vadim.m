function [p c1] = validatepower_Vadim(powertype,saveflag,savetopath,plotmatrix,cir,modifyflag,p1)
%to calculate and graph power comparison
% inputs:
%			feederid : 1 to 5 corresponding to {'355','480','520','909','971'}
%			powertype: 'kvar', 'kw', 'kva' or 'pf' (power factor),'v' (voltage), 'I' (current)
%			saveflag: whether to save figures in /fig folder
%           savetopath
%           plotmatrix: [1 1 1 1 1] to decide what figure(s) to plot. There are 5 of them in total. Fig 1. 2D bar plot. Fig 2. 3D  plot. Fig 3. 3D diff plot. Fig 4: 3D diff plot in per unit. Fig 5: circuit visualizer plot
%           cir : circuit structure for performance speed up. If this is provided, the calculation for it is ignored for performance boost. 
%           glc : similar to cir variable but for given linecode infor
%           modifyflag: whether to use modified version of the circuit
% usage: 
%           validatepower(3,'kvar',0,'',[1 0 0 0 0],c3,glc,1); %most complicated set up 
%           validatepower(3,'kvar'); %most simple set up

if ~exist('saveflag','var')
	saveflag = false;
end

if ~exist('savetopath','var')
	savetopath = 'fig';
end

if ~exist('plotmatrix','var')
	plotmatrix = [1 1 1 1 1];
end

if ~exist('modifyflag','var')
    modifyflag = true;
end

%% create output dir
if ~exist(savetopath,'dir') && ~isempty(savetopath)
	mkdir(savetopath);
end
%% config

% type: 'kvar', 'kw', 'kva' or 'pf' (power factor)
type = powertype;

if strcmp(cir.circuit.Name, 'CENTAUR_12KV')
    fns2 = 'C:\Work\Projects\2013\1665-CEC_Forecasting\Simulations\OpenDSS\Centaur-Jens/Centaur_LoadFlow_130508_Dec2012_Winter_v3.xlsx';
else
    fns2 = 'C:\Work\Projects\2013\1665-CEC_Forecasting\Simulations\OpenDSS\Durox-Jens/Durox_LoadFlow_130730_Dec2012.xlsx';
end

if exist('cir','var')
   c1 = cir;
else
    %% Load power flow data given
    % load data
    d1 = mdb2obj(fns{idx});
    c1 = dssconversion( d1 , glc );
    dsswrite(c1,str{idx},1,str{idx});
end

% utility data
vg = excel2obj(fns2);
vg = vg.('SS_Custom_Load_Flow');
% convert section names to the format we're using for the DSS objects
vg = structconv(vg);
vg.Section_Id_ = regexprep(vg.Section_Id_,' ','_');
vg.Section_Id_ = regexprep(vg.Section_Id_,'\$','_');
vg.Section_Id_ = regexprep(vg.Section_Id_,'\-','_');
vg = structconv(vg);
SECTS = {vg.Section_Id_};
UNI = unique(SECTS);
[aa,bb] = ismember(UNI,{vg.Section_Id_});
vg = vg(bb);

% reorder the sections of utility data
[vgid a] = ismember({vg.Section_Id_},{c1.line.Name});
if(~all(vgid))
	warning('faults:missing',['Some sectionIds are present in the fault data that are missing in the network data:' sprintf('\n\t%s',vg(~a).Section_Id_)]);
	a = a(vgid);
end
all(strcmp({vg(vgid).Section_Id_}',c1.line(a).Name))
vg(~vgid) = '';

% % % % % l = regexprep(c1.line(a).bus2,'(\.\d+)+$','');
% % % % % [vg(vgid).BusID] = deal(l{:});
% % % % % l = regexp({vg(~vgid).Section_Id_}','[^_]+$','match','once');
% % % % % [vg(~vgid).BusID] = deal(l{:});

% get simulated power flow in opendss
t = dssget(p1);
% t = dssget(c1);
% t.Text.Command = 'Set mode=snapshot';
% t.Text.Command = 'Solve';
t.Text.Command = 'Export powers kva power.csv';
pow = excel2obj('power.csv');
pow = pow.power;
pid = ~cellfun(@isempty,regexp({pow.Element}','Line\..+')) & ([pow.Terminal]==1)';
pow = pow(pid);
% clean up power array
psect = strtrim(regexprep({pow.Element},'^(.+\.)',''))';
%pbus = regexprep(psect,'^([a-zA-Z]+\_)','');
%pbus = [regexprep(pbus,'(\_.+)$','')' regexprep(pbus,'^(.+\_)','')'];

% reorder the sections of utility data
[vgid a] = ismember(psect,{c1.line.Name});
if(~all(vgid))
	warning('faults:missing',['Some sectionIds are present in the result data that are missing in the utility data:' sprintf('\n\t%s',psect(~a))]);
	a = a(vgid);
end
all(strcmp(psect(vgid),c1.line(a).Name))
OpenDSS = pow(vgid);
psect = psect(vgid);

% % % % l = regexprep(c1.line(a).bus2,'(\.\d+)+$','');
% % % % [pow(vgid).BusID] = deal(l{:});
% % % % l = regexp(psect(~vgid),'[^_]+$','match','once');
% % % % [pow(~vgid).BusID] = deal(l{:});

% match to buses (simulated data to utility data)
% remove the phase specification from bus name (e.g. bus1.1.2.3 to bus1)
% [vgid pid] = ismember({vg.BusID},{pow.BusID});
[vgid pid] = ismember({vg.Section_Id_},psect);
if(~all(vgid))
	warning('faults:missing','Some busIds are present in lines for the fault data, but not in OpenDSS''s fault output');
	pid = pid(vgid);
end

% % % % p = [vertcat(vg(vgid).Total_Thru_Power_kW_) [pow(pid).P_kW_]'*sqrt(1)];
% % % % q = [vertcat(vg(vgid).Total_Thru_Power_kVAR_) [pow(pid).Q_kvar_]'*sqrt(1)];
p = [vertcat(vg(vgid).Total_Thru_Power_kW_) [OpenDSS(pid).P_kW_]'*sqrt(1)];
q = [vertcat(vg(vgid).Total_Thru_Power_kVAR_) [OpenDSS(pid).Q_kvar_]'*sqrt(1)];

LinesUtility = {vg(vgid).Section_Id_}';
LinesOpenDSS = {OpenDSS(pid).Element}';

% recalculate p based on config param
switch lower(type)
	case 'V'
		p = [vertcat(vg(vgid).Volts_Out) [pow(pid).P_kW_]'*sqrt(3)];
	case 'I'
		
	case 'kw'
	case 'kvar'
		p = q;
	case 'kva'
		p = (p.^2 + q.^2).^.5;
	case 'pf'
		p = p./(p.^2 + q.^2).^.5;
end
% distance from feeder
dist = vertcat(vg(vgid).Total_distance_ft)/1000;%in kft
% lookup xy coordinates
% % % % [xy xy] = ismember({vg(vgid).BusID}',c1.buslist.id);
% % % % xy = c1.buslist.coord(xy,:);
% and get a sort order (now named b) for sorting by distance
[order order] = sort(dist);
TEST = {};

%% plot real power with difference bar
if plotmatrix(1)
	figure; hold on;
	for i=1:length(dist)
		if(isnan(p(order(i),1))), continue; end
		%plot([dist(order(i)) dist(order(i))],abs(powerfactor(order(i),:)),'r');
		plot([dist(order(i)) dist(order(i))],p(order(i),:),'r');
        test = {LinesUtility(order(i),:), LinesOpenDSS(order(i),:),p(order(i),1),p(order(i),2)};
        TEST = [test;TEST];
	end
	%h = plot(dist(order),abs(powerfactor(order,1)),'x',dist(order),abs(powerfactor(order,2)),'.');
	h = plot(dist(order),p(order,1),'x',dist(order),p(order,2),'.');
	legend(h,{'Source Data','OpenDSS'});
	xlabel('Distance, kft')
	ylabel(type);
    grid;
	set(gcf,'nextplot','new')
	if saveflag 
		saveas(gcf,[savetopath '/feeder' str{idx} '_2D_P_' type '.fig']);
	end
end

% plot real power difference between utility data and opendss data
% figure; hold on;
% h = plot(dist(order),p(order,1)-p(order,2),'x',dist(order),p1(order,1)-p1(order,2),'.');
% legend(h,{'Source Data - OpenDSS Data'});
% xlabel('Distance, kft')
% ylabel('Difference in Real Power, kw');
% set(gcf,'nextplot','new')

% difference bars in 3D
% [xy xy] = ismember({vg(vgid).BusID},c1.buslist.id);
% xy = c1.buslist.coord(xy,:);
% figure; hold on;
% for i=1:length(dist)
% 	if(isnan(p(order(i),1))), continue; end
% 	%plot([dist(order(i)) dist(order(i))],abs(powerfactor(order(i),:)),'r');
% 	plot3([xy(i,1) xy(i,1)],[xy(i,2),xy(i,2)],abs(p(order(i),:)),'r');
% end
% %h = plot(dist(order),abs(powerfactor(order,1)),'x',dist(order),abs(powerfactor(order,2)),'.');
% h = plot3(xy(:,1),xy(:,2),abs(p(order,1)),'x',xy(:,1),xy(:,2),abs(p(order,2)),'.');
% legend(h,{'Source Data','OpenDSS'});
% xlabel('Distance, kft')
% zlabel('Real Power, kw');
% set(gca,'looseinset',get(gca,'tightinset'));
% set(gcf,'nextplot','new')

% real power: plot per unit comparison 3D
% Draw the circuit in 3D with difference on the z-axis
% lookup lines and map the data onto them
% m = regexprep({c1.line.bus1}','(\.\d+)+$','');
% n = regexprep({c1.line.bus2}','(\.\d+)+$','');
% [ism_ locm] = ismember(m,{pow(pid).BusID});
% [isn_ locn] = ismember(n,{pow(pid).BusID});
% mask = ism_ & isn_;
% locm = locm(mask);
% locn = locn(mask);
% % collect the data for plotting
% % (x1,y1,z1,x2,y2,z2,nan,nan,nan) for each line, then reshape to be able to
% % plot as a single object
% %dat_ = [xy(locm,:),pow(locm).P_kW_,xy(locn,:),cdiff(locn),nan(length(locm),3)];
% %dat_ = [xy(locm,:),[pow(locm).P_kW_]',xy(locn,:),[pow(locn).P_kW_]',nan(length(locm),3)];
% %pow = [vg(vgid).Into_c_kVA]';
% 
% %%
% if plotmatrix(2)
% 	dat_ = [xy(locm,:),p(locm,2),xy(locn,:),p(locn,2),nan(length(locm),3)];
% 	dat_2 = [xy(locm,:),p(locm,1),xy(locn,:),p(locn,1),nan(length(locm),3)];
% 	dat_ = reshape(dat_',3,numel(dat_)/3)';
% 	dat_2 = reshape(dat_2',3,numel(dat_2)/3)';
% 	figure; hold on
% 	plot3(dat_(:,1),dat_(:,2),dat_(:,3))
% 	plot3(dat_2(:,1),dat_2(:,2),dat_2(:,3),'r')
% 	xlabel('X'); ylabel('Y'); zlabel(type);
% 	legend({'OpenDSS','Source Data'});
% 	set(gca,'looseinset',get(gca,'tightinset'));
% 	set(gcf,'nextplot','new')
% 	if saveflag
% 		saveas(gcf,[savetopath '/feeder' str{idx} '_3D_P_' type '.fig']);
% 	end
% end
% 
% % 3D plot kw diff
% if plotmatrix(3)
% 	d_ = -p(:,1)+p(:,2);
% 	dat_ = [xy(locm,:),d_(locm),xy(locn,:),d_(locn),nan(length(locm),3)];
% 	dat_ = reshape(dat_',3,numel(dat_)/3)';
% 	figure; hold on
% 	plot3(dat_(:,1),dat_(:,2),dat_(:,3))
% 	xlabel('X'); ylabel('Y'); zlabel(['\Delta P, ' type]);
% 	legend({'OpenDSS - Source Data'});
% 	set(gca,'looseinset',get(gca,'tightinset'));
% 	set(gcf,'nextplot','new')
% 	if saveflag
% 		saveas(gcf,[savetopath '/feeder' str{idx} '_3D_dP_' type '.fig']);
% 	end
% end
% % 3D plot pu diff
% if plotmatrix(4)
% 	d_ = (p(:,2)-p(:,1))./p(:,1);
% 	dat_ = [xy(locm,:),d_(locm),xy(locn,:),d_(locn),nan(length(locm),3)];
% 	dat_ = reshape(dat_',3,numel(dat_)/3)';
% 	figure; hold on
% 	plot3(dat_(:,1),dat_(:,2),dat_(:,3))
% 	xlabel('X'); ylabel('Y'); zlabel(['\Delta P(' type '), pu']);
% 	legend({'(OpenDSS - Source Data) in pu'});
% 	set(gca,'looseinset',get(gca,'tightinset'));
% 	set(gcf,'nextplot','new')
% 	if saveflag
% 		saveas(gcf,[savetopath '/feeder' str{idx} '_3D_dP_' type '.fig']);
% 	end
% end
% 
% if plotmatrix(5)
% 	circuitVisualizer(c1)
% 	if saveflag
% 		saveas(gcf,[savetopath '/feeder' str{idx} '_2D_visual.fig']);
% 	end
end