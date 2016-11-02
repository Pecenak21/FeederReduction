% Procedure for running OpenDSS conversion code

flgTest=0; % set to '1' to select some reduced versions of the files that describe the system


%% Get data 
cd('D:\Sync Software\Other Code\SystemConversion\SCE_Centaur');
curdir = pwd;
fns = {[curdir '\input\Centaur12kV_Load.txt'];...
        [curdir '\input\Centaur12kV_Network.txt'];...
        [curdir '\input\Centaur12kV_Equipment.txt']};
flgChopFile=0;
% fns={}; % activate if larger files are already chopped

% load excel data
d1 = cyme2obj_SCE(fns,flgChopFile,flgTest);
% d2 = excel2obj(fns{2});
% d3 = excel2obj(fns{3});
% d4 = excel2obj(fns{4});
% d5 = excel2obj(fns{5});

% given linecode
% glc = excel2obj( 'linecode.xlsx' );
% glc = glc.LineCode;
glc1=d1.LineCode;

%% or just load data from disk (if available)
% load exceldata.mat

%% Feed data in conversion engine
fprintf('converting excel to opendss ...');tic
c1 = dssconversion_CYME(d1,glc1,'355',flgTest);
toc

%% Write out OpenDSS files
fprintf('writing ...');tic
p1 = dsswrite(c1,'355',1,'355');
toc

faultstudy_CYME(c1,d1,p1,0);

%% Read back
fprintf('parsing ...');tic
[tc1 cm1] = dssparse(p1);
toc

%% compare opendss struct results
dssstructcmp(c1,tc1);

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