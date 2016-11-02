%% Output Options
strOutput='dss';
% strOutput='emtp';

%% Input Options

% Select system
% strSystem='HydroOttawa';
% strSystem='SDGE_all'; 
% strSystem='SDGE_520';
strSystem='SCE_Centaur';
% strSystem='SCE_Durox';

% map system to software
cellSource_cyme=[{'HydroOttawa'};{'SCE_Centaur'};{'SCE_Durox'}];
cellSource_synergee=[{'SDGE_all'};{'SDGE_520'}];

% Select format of input system files
if ismember(strSystem,cellSource_cyme)
    strSource='cyme';
elseif ismember(strSystem,cellSource_synergee)
    strSource='synergee';
else
    strSource='cyme'; %default
end


flgTest=0; % set to '1' to select some reduced versions of the files that describe the system
% flgChopFile=0; % chop large CYME files into smaller files

% path to system files and default units
if strcmp(strSystem,'HydroOttawa')
    strDirectory='D:\Sync Software\Other Code\SystemConversion\HydroOttawa';
    strDirectory_output=strDirectory;
    % HydroOttawa info for CYME data: length unit is ft and the impedances are ohms/mile
    NetworkID='';
    Unit_length='ft';
    Unit_impedance='OhmsPerMile';
    Unit_geometry_length='m';
    MVA='100'; % Need to assume 100 MVA as base to match customer-provided fault data
elseif strcmp(strSystem,'SCE_Centaur')
    strDirectory='D:\Sync Software\Other Code\SystemConversion\SCE_Centaur';
    strDirectory_output=strDirectory;
    % From Sunil's email on 5/23/2013:
    % Spacing is meters, radius is centimeters, and resistance is ohm/km
    % Previously, Sunil told us on the phone that line length is in feet.
    % Better double check this information as it seems inconsistent with 
    % the other units, which are metric.
    NetworkID='CENTAUR_12KV';
    Unit_length='ft';
    Unit_impedance='OhmsPerKMeter';
    Unit_geometry_length='m';
    MVA='';
elseif strcmp(strSystem,'SCE_Durox')
    strDirectory='D:\Sync Software\Other Code\SystemConversion\SCE_Centaur'; % Centaur and Durox system provided in one set of files
    strDirectory_output=strDirectory;
    NetworkID='DUROX_12KV';
    Unit_length='kft';
    Unit_impedance='OhmsPerMile';
    Unit_geometry_length='m';
    MVA='';
elseif strcmp(strSystem,'SDGE_all')
    strDirectory='D:\Sync Software\Other Code\SystemConversion\SDGE';
    strDirectory_output=strDirectory;
    fns = {[strDirectory '/input/355ForENERNEX.xlsx'];...
            [strDirectory '/input/480ForENERNEX.xlsx'];...
            [strDirectory '/input/520ForENERNEX.xlsx'];...
            [strDirectory '/input/909ForENERNEX.xlsx'];...
            [strDirectory '/input/971ForENERNEX.xlsx']};
    Unit_length='ft';
    Unit_impedance='OhmsPerKFeet';
    MVA='';
elseif strcmp(strSystem,'SDGE_520')
    strDirectory='D:\Sync Software\Other Code\SystemConversion\SDGE'; 
    strDirectory_output=strDirectory;
    NetworkID='';
    fns = {[strDirectory '\input\520ForENERNEX.xlsx']};
    Unit_length='ft';
    Unit_impedance='OhmsPerKFeet';
    MVA='';
end

% Convert impedances to ohms per unit length
cf=ConvertUnit(Unit_impedance,Unit_length);



    
cd(strDirectory);    
curdir = pwd;





%% Reading Input Files
% specify locations of source files here
% process input data and put into a structure that has the complete system information

if strcmp(strSource,'cyme')
    o=cymeread_dir([strDirectory '\input']);
    if flgTest
        fid = fopen([strDirectory '\input\TestSystem.txt'], 'r');
        y=textscan(fid,'%s');
        fclose(fid);
        o.TestSystem = y{:};
    end
    d1 = cyme2obj(o,flgTest,NetworkID);
    glc=d1.LineCode;
    gwd=d1.WireData;
    glg=d1.LineGeometry;
    % Add system-specific defaults (information that is not contained in the source files
    [glc.ConversionFactor]=deal(cf);
    [glc.Unit_length]=deal(Unit_length);
    [glc.Unit_impedance]=deal(Unit_impedance);
    [glg.Unit_length]=deal(Unit_geometry_length);
%     if ~strcmp(MVA,'')
        n_sources=length(d1.Sources);
        for i_=1:n_sources
            d1.Sources(i_).MVA=MVA;
        end
%     end
elseif strcmp(strSource,'synergee') 
    if strcmp(strSystem,'SDGE_all')
        o1 = excel2obj(fns{1});
        o2 = excel2obj(fns{2});
        o3 = excel2obj(fns{3});
        o4 = excel2obj(fns{4});
        o5 = excel2obj(fns{5});
        d1=synergee2cyme(o1);
        d2=synergee2cyme(o2);
        d3=synergee2cyme(o3);
        d4=synergee2cyme(o4);
        d5=synergee2cyme(o5);
    elseif strcmp(strSystem,'SDGE_520')
        % d1 = mdb2obj(fns{1}); 
            % This function load the data directly from the access db,
            % which has the advantage that the data does not have to be
            % exported to Excel. However, this only works with the 32 bit
            % version of Matlab (not the 64 bit version), which is a bug.
        o = excel2obj(fns{1});
        glc = excel2obj( 'linecode.xlsx' );
        o.LineCode=glc.LineCode;
        ob=synergee2cyme(o);
        d1 = cyme2obj(ob,flgTest,NetworkID);
        %     if ~strcmp(MVA,'')
        n_sources=length(d1.Sources);
        for i_=1:n_sources
            d1.Sources(i_).MVA=MVA;
        end
        %     end
    end
    
    glc = d1.LineCode;
    % Add system-specific defaults (information that is not contained in the source files
    [glc.ConversionFactor]=deal(cf);
    [glc.Unit_length]=deal(Unit_length);
    [glc.Unit_impedance]=deal(Unit_impedance);
    glc=structconv(glc);
end


% keyboard



%% Feed data in conversion engine
if strcmp(strSource,'cyme') && strcmp(strOutput,'dss')
    fprintf('converting to Matlab structure in opendss format ...');tic
    c1 = dssconversion_CYME( d1 , glc , gwd, glc);
    toc
elseif strcmp(strSource,'synergee') && strcmp(strOutput,'dss')
    if strcmp(strSystem,'SDGE_all')
        fprintf('converting to Matlab structures in opendss format ...');tic
        c1 = dssconversion( d1 , glc );
        c2 = dssconversion( d2 , glc );
        c3 = dssconversion( d3 , glc );
        c4 = dssconversion( d4 , glc );
        c5 = dssconversion( d5 , glc );
        toc
    else 
        fprintf('converting to Matlab structure in opendss format ...');tic
%         c1 = dssconversion( d1 , glc );
        c1 = dssconversion_CYME( d1 , glc );
        toc
    end
elseif strcmp(strSource,'cyme') && strcmp(strOutput,'emtp')
    
elseif strcmp(strSource,'synergee') && strcmp(strOutput,'emtp')
        
end

% keyboard

%% Tweak circuit
if strcmp(strSystem,'SDGE_520')
    c1 = circuit_changes_SDGE_520(c1);
end


%% Write out files

if strcmp(strOutput,'dss')
    cd(strDirectory_output);    
    curdir = pwd;

    if strcmp(strSystem,'SDGE_all')
        fprintf('writing ...');tic
        p1 = dsswrite(c1,'355',1,'355');
        p2 = dsswrite(c2,'480',1,'480');
        p3 = dsswrite(c3,'520',1,'520');
        p4 = dsswrite(c4,'909',1,'909');
        p5 = dsswrite(c5,'971',1,'971');
        toc
    else 
        fprintf('writing ...');tic
        p1 = dsswrite(c1,strSystem,1,strSystem);
        toc
    end
elseif strcmp(strOutput,'emtp')
    
end

keyboard



%% Do other stuff
% at this point the system should be converted
% write here what you want to do with the converted system


% faultstudy_CYME(c1,d1,p1,0);
% runDSS;

% 
% %% Read back
% fprintf('parsing ...');tic
% [tc1 cm1] = dssparse(p1);
% [tc2 cm2] = dssparse(p2);
% [tc3 cm3] = dssparse(p3);
% [tc4 cm4] = dssparse(p4);
% [tc5 cm5] = dssparse(p5);
% toc
% 
% %% compare opendss struct results
% dssstructcmp(c1,tc1);
% dssstructcmp(c2,tc2);
% dssstructcmp(c3,tc3);
% dssstructcmp(c4,tc4);
% dssstructcmp(c5,tc5);
% 
% %% run simulation in matlab and get properties of a load out (for example)
% [l1 circuit]= dssget(c1,c1.load(1));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% TEST with feeder 13 
% % cmd = 'show voltages LN nodes';
% 
% % run simulation on original file provided in OpenDSS folder
% filename = 'IEEE13Nodeckt.dss';
% [a1 b1] = dssget(filename);
% 
% % convert to opendss struct and run it
% [tc cmds] = dssparse(filename);
% p2f = dsswrite(tc,'ieee13',0,'ieee13',cmds);
% 
% [a2 b2] = dssget(p2f);
% dssstructcmp(b1,b2);
% 
% %% test IEEE bus30 circuit
% % run simulation on original file provided in OpenDSS folder
% filename = 'IEEETestCases\IEEE 30 Bus\Run_IEEE30.DSS';
% [a1 b1] = dssget(filename);
% 
% % convert to opendss struct and run it
% [tc cmds] = dssparse(filename);
% p2f = dsswrite(tc,'ieee30',0,'ieee30',cmds);
% 
% [a2 b2] = dssget(p2f);
% dssstructcmp(b1,b2);
% 
% %% These should be working but you might consider if you really want to run
% % them. They take a while.
% % Test EPRITestCase ckt24
% [c24 cmd24] = dssparse('EPRITestCircuits\ckt24/Run_Ckt24.dss');
% p2f = dsswrite(c24,'ckt24',1,'ckt24',cmd24);
% 
% dssget('EPRITestCircuits\ckt24/Run_Ckt24.dss');
% dssget(p2f);
% 
% % Test EPRITestCase ckt7
% [c7 cmd7] = dssparse('EPRITestCircuits\ckt7/RunDSS_Ckt7.dss');
% p2f = dsswrite(c7,'ckt7',1,'ckt7',cmd7);
% 
% dssget('EPRITestCircuits\ckt7/RunDSS_Ckt7.dss');
% dssget(p2f);
% 
% % Test EPRITestCase ckt5
% [c5 cmd5] = dssparse('EPRITestCircuits\ckt5/Run_Ckt5.dss');
% p2f = dsswrite(c5,'ckt5',1,'ckt5',cmd5);
% 
% dssget('EPRITestCircuits\ckt5/Run_Ckt5.dss');
% dssget(p2f);
% 
% %% TODO: does not work on following circuits yet
% %% Test 8500 node feeder
% [c85 cmd85] = dssparse('IEEETestCases/8500-Node/Run_8500Node.dss');
% 
% %% No need to test following parts
% %% TEST: Find unknown linecodes
% % Line resistance unit: ohm/1000ft
% dat = {d1, d2, d3, d4, d5};
% unknownCond = [];
% for k = 1:length(dat)
%     d = dat{k};
%     lc = excel2obj( [curdir '\linecode.xlsx'] );
% 
%     % Check for missing linecode
%     cond = unique([d.InstSection.PhaseConductorId; d.InstSection.NeutralConductorId]);
%     fcond = lc.LineCode.ConductorId;
% 
%     id = [];
%     for i = 1:length(cond)
%         if sum( strcmp(cond{i}, fcond) ) < 1
%             id = [id i];
%         end
%     end
%     unknownCond = unique([unknownCond; cond(id)]);
% end
% 
% %% TEST: Draw Loads' position
% close all
% drawSections('generate',d);
% drawSections(d.InstSection.SectionId,'b');
% drawSections(d.Loads.SectionId,'r',1);
% hold on;
% 
% % node indexing based on names
% % build a hash lookup table for node coords
% 	for i=1:length(d.Node.NodeId)
% 		nodeXY.(['n' d.Node.NodeId{i}]) = [d.Node.X{i} d.Node.Y{i}];
% 	end
% 
% % section indexing based on names
% % build a hash lookup table for sections
% for i=1:length(d.InstSection.SectionId)
%     id = ['s' d.InstSection.SectionId{i}];
%     id(id==' ') = [];
%     sectionH.(id) = struct('FromNodeId',d.InstSection.FromNodeId{i},'ToNodeId',d.InstSection.ToNodeId{i});
% end
% %
% for i=1:length(d.Loads.SectionId)
% 	id = ['s' d.Loads.SectionId{i}];
% 	id(id==' ') = [];
% 	section = sectionH.(id);
% 	coords = vertcat(nodeXY.(['n' section.ToNodeId]));
% 	plot(coords(1),coords(2),'o');
% end