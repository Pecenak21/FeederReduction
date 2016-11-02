function output = dssconversion_CYME(circuitdata, linecodedata, gwd, circuitname)
% This conversion engine works for data from CYME only. 
% Modification is needed if other source of data is used.
% Inputs:
%			circuitdata: structure with CYME data
%			linecodedata: linecode data or path to linecode excel file
%			wiredata: used for lines specified by line geometry
%			circuitname: (optional). Default: 'newcircuit'

% Process inputs
if isstruct(circuitdata)
	d = circuitdata;
else
	error('Invalid circuitdata type');
end

if isstruct(linecodedata)
	glc = linecodedata;
elseif isa(linecodedata, 'char')
	fprintf('Reading linecodes... '); t_ = tic;
	glc = excel2obj(linecodedata);
	glc = glc.LineCode;
	toc(t_);
else
	error('Invalid linecodedata type');
end


% phasemap.(d.Section(1).Phase1) = '.1';
% phasemap.(d.Section(1).Phase2) = '.2';
% phasemap.(d.Section(1).Phase3) = '.3';
% phasemap.N = '.0';
phasemap.X = '.1';
phasemap.Y = '.2';
phasemap.Z = '.3';
phasemap.N = '.0';

% keep track of fieldnames we handle
% when we're all done we'll report unhandled fields to the user
fns_handled = {};

%% Convert linecode to OpenDSS linecode
fprintf('Converting linecodes... '); t_ = tic;
lc=Convert2DSS_LineCode(glc);
toc(t_);

%% Convert wiredata to OpenDSS wiredata
if exist('gwd')
    fprintf('Converting wiredata... '); t_ = tic;
    wd=Convert2DSS_WireData(gwd);
    toc(t_);
end


% %% Convert line geometry to OpenDSS line geometry
% fprintf('Converting line geometry... '); t_ = tic;
% lg=Convert2DSS_LineGeometry(glg);
% toc(t_);

%% Convert sections to OpenDSS objects
fprintf('Converting sections... '); t_ = tic;
% UG and OH lines
% [flg,index]=ismember(d.Section.SectionType,{'UG'});
i_line=1;
i_geometry=1;
% i_load=1;
i_capacitor=1;
i_transformer=1;
lg=dsslinegeometry;
lg_name={''};
for i=1:length(d.Section)
    % go through each section, look at the section type, and create the
    % appropriate project (line specified by geometry, capacitor, or
    % transformer
    if strcmp(d.Section(i).SectionType,'OH') || strcmp(d.Section(i).SectionType,'UG') || strcmp(d.Section(i).SectionType,'geometry')
        if strcmp(d.Section(i).SectionType,'geometry')
            LineGeometry={[d.Section(i).SpacingID '__' d.Section(i).CondID_A '__' d.Section(i).CondID_B '_' d.Section(i).CondID_C '__' d.Section(i).CondID_N]};  
            if ~ismember(LineGeometry,lg_name) % add to geometry library if this line does not already exist
                LineGeometry_conv=structconv(d.LineGeometry);
                SpacingID=LineGeometry_conv.ID;
                [flg index]=ismember({d.Section(i).SpacingID},SpacingID);
                lg(i_geometry)=Convert2DSS_LineGeometry(d.Section(i),d.LineGeometry(index),LineGeometry);
                i_geometry=i_geometry+1;
                lg_name=[lg_name LineGeometry];
            end
        end
        l(i_line)=Convert2DSS_Lines(d.Section(i),lc);
        i_line=i_line+1;
%     elseif strcmp(d.Section(i).SectionType,'geometry')
%         LineGeometry=[d.Section(i).LineGeometry];
%         if ~ismember({LineGeometry},lg.Name)
%             lg(i_geometry)=Convert2DSS_Geometry(d.Section(i),wd,lg);
%             i_geometry=i_geometry+1;
%         end
    elseif strcmp(d.Section(i).SectionType,'switch')
        l(i_line)=Convert2DSS_Lines(d.Section(i),lc);
        i_line=i_line+1;
    elseif strcmp(d.Section(i).SectionType,'other') % model other sections as switches
        l(i_line)=Convert2DSS_Lines(d.Section(i),lc);
        i_line=i_line+1;
    elseif strcmp(d.Section(i).SectionType,'capacitor')
        % JS: Not sure how phase connection of caps is being determined,
        % need to revisit, probably need to fold in .network.section info
        [cp(i_capacitor) capcon(i_capacitor)]=Convert2DSS_Capacitors(d.Section(i));
        i_capacitor=i_capacitor+1;
        d_Section=structconv(d.Section);
        
        % Caps in the CYME system have a unique section ID, i.e., the caps
        % are not connected to the system. Caps in the SynerGEE system have
        % a section ID that is also a section ID of a line object, i.e.,
        % the caps are connected to the sytem via that line object. 
        % For the CYME system, we need to connect the caps to the system
        % via a line. We do this here by inserting a 'switch' line (OpenDSS
        % 'switch' property is set to 'yes').
        i_=find(ismember(d_Section.SectionID,d.Section(i).SectionID));
        if length(i_)==1 % The assumption is that if the capacitor Section ID occurs only once, then we need to create a switch line that connect the capacitor to the system  
            l(i_line)=Convert2DSS_Lines(d.Section(i),lc);
            i_line=i_line+1;
        end        
    elseif strcmp(d.Section(i).SectionType,'transformer')
        tr(i_transformer)=Convert2DSS_Transformers(d.Section(i));
        i_transformer=i_transformer+1;
    end
end



% l=Convert2DSS_Lines(d.Section,lc);
fns_handled{end+1} = 'Section';
if i_line>1
    fns_handled{end+1} = 'Line';
end
if i_geometry>1
    fns_handled{end+1} = 'Geometry';
end
% if i_load>1
%     fns_handled{end+1} = 'Load';
% end
if i_capacitor>1
    fns_handled{end+1} = 'Capacitor';
end
if i_transformer>1
    fns_handled{end+1} = 'Transformer';
end
toc(t_);


%% Convert loads to OpenDSS objects

if isfield(d,'Loads')
    fprintf('Converting loads... '); t_ = tic;
    i_load=1;
    i_ge=1;
    for i=1:length(d.Loads)
        [ld_temp ge_temp]=Convert2DSS_Loads(d.Loads(i));
        for i_=1:length(ld_temp)
            ld(i_load)=ld_temp(i_);
            i_load=i_load+1;
        end
        for i_=1:length(ge_temp)
            ge(i_ge)=ge_temp(i_);
            i_ge=i_ge+1;
        end
    end
    fns_handled{end+1} = 'Load';
    toc(t_);
end



%% Switches
% Benandu's code for SynerGEE handle switches separately and create a
% 'switch' OpenDSS file. I use a different approach by incorporating
% switches as lines with the 'switch' parameter set to 'yes' if the switch
% is closed. Tested with CYME, needs to be tested with SynerGEE.

% if(isfield(d,'Switches'))
% 	fprintf('Converting switches... '); t_ = tic;
% 	sws = d.Switches;
% 	if(any([sws.Automated]))
%         warning('dssconversion:Unimplemented','Automated switching is currently not implemented; if this feature is important to your project, please add that to the code!');
%     end
% %     if(any([sws.SwitchIsTie]))
% %         warning('dssconversion:SwitchIsTie:NeedToCheck','Implemented SwitchIsTie as "lock" parameter of switch (best guess). You may need to double check and see if this is what you want in your circuit!');
% %     end
%     % Determine which unique ids are actually unique (used for naming)
%     % first get them all
%     useuid = {sws.DeviceNumber};
%     % then remove one instance of each
%     [tmp i] = unique(useuid);
%     useuid(i) = [];
%     % and then don't use any of the remaining items (i.e. anything that was
%     % present more than once)
%     useuid = ~ismember({sws.DeviceNumber},useuid);
%     for i = 1:length(sws)
%         sw(i)=Convert2DSS_Switches(sws(i),useuid(i));    
%     end
% 	toc(t_);
% 	fns_handled{end+1} = 'Switches';
% end



%% Reclosers
if(isfield(d,'Reclosers'))
	fprintf('Converting reclosers ... '); t_ = tic;
	rcs = d.Reclosers;
    for i = 1:length(rcs)
        rc(i)=Convert2DSS_Reclosers(rcs(i));    
    end
	toc(t_);
	fns_handled{end+1} = 'Reclosers';
end

%% convert fuses
if(isfield(d,'Fuses'));
	fprintf('Converting fuses... '); t_ = tic;
	fs = d.Fuses;
    if ~isempty(fs)
        for i = 1:length(fs)
            f(i)=Convert2DSS_Fuses(fs(i));
        end
    end
	toc(t_);
	fns_handled{end+1} = 'Fuses';
end

%% convert capacitors
if(isfield(d,'Capacitors'));
	fprintf('Converting capacitors... '); t_ = tic;
	cps = d.Capacitors;
    if ~isempty(cps)
        for i = 1:length(cps)
            [cp(i) capcon(i)]=Convert2DSS_Capacitors(cps(i));
        end
    end
	toc(t_);
	fns_handled{end+1} = 'Capacitors';
end

%% Transformers
% transformers are convered in 'Convert sections to OpenDSS objects'
% here we only define trs, which is necessary for dealing with regulators
% in the next section

% js: not sure how the zero sequence impedances plays in, I suppose Z1 and
% Z0 are assumed to be equal
%
% js: do the X/R and X0/R0 ratios matter for OpenDSS?
% 
% js: apparently, OpenDSS allows to specify the properties of each winding
% separatetly or, alternatively, the properties of all windings via arrays,
% below the properties of the windings are specified as arrays, but the
% default values for the individual windings are still retained, which
% might confuse OpenDSS because values are not consistent (certainly 
% confuses the user), better not use default values or set them to zero if 
% arrays are used

if(isfield(d,'Transformers'))
% 	fprintf('Converting transformers... '); 
%     t_ = tic;
	trs=d.Transformers;
% 	for i = 1:length(trs)
%         tr(i)=Convert2DSS_Transformers(trs(i));
%     end
%     %remove lines that represent these transformer(s)
%     [val l_idx] = ismember(dataclean({d.Transformers.SectionID},'name'),{l.Name});
% 	l(nonzeros(l_idx)) = [];
%     toc(t_);
% 	fns_handled{end+1} = 'Transformers';
else
    trs=[];
end




%% Only translated for SynerGEE. Need to modify to be able to convert regulators from CYME.
% %% Regulator Controllers (should be defined after transformers)
% % ASSUMPTION: All raw transformer raw data is stored in 'trs' variable
if(isfield(d,'Regulators'))
	fprintf('Converting regulator controllers... '); t_ = tic;
    rs=d.Regulators;
    for i = 1:length(rs)
        [r(i) tr(i)]=Convert2DSS_Regulators(rs(i),trs,phasemap);
    end
    %remove lines that represent these regulator(s)
    [val l_idx] = ismember(dataclean({d.Regulators.SectionId},'name'),{l.Name});
	l(nonzeros(l_idx)) = [];
    toc(t_);
    fns_handled{end+1} = 'Regulators';
end 
    
  



%% Circuit/Feeder
% js: potential issue are
% - is more than one feeder source allowed if treated as substation?
% - how about .bus2 parameter? Assumed to be grounded if not specified?
% - OpenDSS gets confused if MVA,X1R1 is given (default values) and, the
%   alternative input method, R1,X1,R0,X0 is used, better not specify default
%   values
% - OperatingVoltage1, OperatingVoltage2, OperatingVoltage3 form CYME data
%   not used, perhaps .pu = OperatingVoltage/.basekv (assuming balance)

% % According to the SynerGEE documentation, there are two kinds of sources
% % in synergee.  A "feeder" is used when substation data is unavailable, and
% % is "usually setup to represent some point near the secondary of the
% % substation transformer."  A "substation" had historically included the
% % transformer as well, but that data has been moved to the primary
% % transfomers table when applicable, leaving the "substation" very similar
% % to the feeder.  Not sure how "substation" is stored in a table since we
% % don't have an example of one.


fprintf('Converting sources... '); t_ = tic;
fi_source=structconv(d.Sources);
fi = generatesourcedata(fi_source);
n_circuit=1;  % OpenDSS only allows one circuit-defining source, source #1 is Centaur, source #2 is Durox

% generate one circuit-defining source and 0 to many voltage sources
% circuit with multiple sources not tested, yet
[cir vs]=Convert2DSS_Sources(fi,n_circuit);
toc(t_);
fns_handled{end+1} = 'Feeder';


%% Bus list with locations
if(isfield(d,'Nodes'))
	fprintf('Generating buslist... '); t_ = tic;

	output.buslist.id = dataclean({d.Nodes.NodeID},'name')';
    if strcmp(class(d.Nodes(1).CoordX),'double')
        output.buslist.coord = [vertcat(d.Nodes.CoordX), vertcat(d.Nodes.CoordY)];
    else
        output.buslist.coord = str2double([{d.Nodes.CoordX}', {d.Nodes.CoordY}']);
    end

	fns_handled{end+1} = 'Nodes';

	toc(t_);
end

%% Check what fields are left
fns_handled = [fns_handled {'InstMymLargeCust','InstMymLoads','InstMymMeter','InstViews','SAI_Control'}]; % to specify others all at once
fn_unknown = setdiff(fieldnames(d),fns_handled);
disp('Unknown object types:');
disp(fn_unknown);

%% Output
voltages = [];
if(exist('cir','var'))
	output.circuit = cir;
	voltages = [voltages cir.basekv];
end
if(exist('capcon','var'))
	output.capcontrol = capcon;
end
if(exist('ld','var'))
	output.load = ld;
	voltages = [voltages ld.Kv];
end
if(exist('ge','var'))
	output.generator = ge;
end
if(exist('l','var'))
	output.line = l;
end
if(exist('lc','var'))
	output.linecode = lc;
end
if(exist('sw','var'))
	output.switch = sw;
end
if(exist('r','var'))
	output.regcontrol = r;
end
if(exist('f','var'))
	output.fuse = f;
end
if(exist('rc','var'))
	output.recloser = rc;
end
if(exist('cp','var'))
	output.capacitor = cp;
	voltages = [voltages cp.Kv];
end
if(exist('tr','var')) && ~isempty(tr)
	output.transformer = tr;
	voltages = [voltages tr.kVs];
end
if(~isempty(voltages))
	output.basevoltages = sort(unique(voltages),'descend');
end
