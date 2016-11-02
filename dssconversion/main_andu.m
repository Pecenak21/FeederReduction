%% Testing Script
% Get data 
curdir = [pwd '/dssconversion'];
% fns = strcat('custdata/', ...
% 		{'355','480','520','909','971'}','ForENERNEX.mdb');
% 
% % load excel data
% d1 = mdb2obj(fns{1});
% d2 = mdb2obj(fns{2});
% d3 = mdb2obj(fns{3});
% d4 = mdb2obj(fns{4});
% d5 = mdb2obj(fns{5});

fns = {[curdir '/custdata/355ForENERNEX.xlsx'];...
        [curdir '/custdata/480ForENERNEX.xlsx'];...
        [curdir '/custdata/520ForENERNEX.xlsx'];...
        [curdir '/custdata/909ForENERNEX.xlsx'];...
        [curdir '/custdata/971ForENERNEX.xlsx']};

% load excel data
d1 = excel2obj(fns{1});
d2 = excel2obj(fns{2});
d3 = excel2obj(fns{3});
d4 = excel2obj(fns{4});
d5 = excel2obj(fns{5});

% given linecode
glc = excel2obj( [curdir 'dssconversion\linecode.xlsx'] );
glc = glc.LineCode;

% or just load data from disk (if available)
% load exceldata.mat

%% Feed data in conversion engine
fprintf('converting excel to opendss ...');tic
c1 = dssconversion( d1 , glc );
c2 = dssconversion( d2 , glc );
c3 = dssconversion( d3 , glc );
c4 = dssconversion( d4 , glc );
c5 = dssconversion( d5 , glc );
toc

%% Write out OpenDSS files
p = 'toSDGE';
fprintf('writing ...');tic
p1 = dsswrite(c1,'355',1,[p '/355']);
p2 = dsswrite(c2,'480',1,[p '/480']);
p3 = dsswrite(c3,'520',1,[p '/520']);
p4 = dsswrite(c4,'909',1,[p '/909']);
p5 = dsswrite(c5,'971',1,[p '/971']);
toc

%% Read back
fprintf('parsing ...');tic
[tc1 cm1] = dssparse(p1);
[tc2 cm2] = dssparse(p2);
[tc3 cm3] = dssparse(p3);
[tc4 cm4] = dssparse(p4);
[tc5 cm5] = dssparse(p5);
toc

%% compare opendss struct results
dssstructcmp(c1,tc1);
dssstructcmp(c2,tc2);
dssstructcmp(c3,tc3);
dssstructcmp(c4,tc4);
dssstructcmp(c5,tc5);

%% run simulation in matlab and get properties of a load out (for example)
[l1 circuit]= dssget(c1,c1.load(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST with feeder 13 
% cmd = 'show voltages LN nodes';

% run simulation on original file provided in OpenDSS folder
filename = 'dssconversion\13Bus\IEEE13Nodeckt.dss';
[a1 b1] = dssget(filename);

% convert to opendss struct and run it
[tc cmds] = dssparse(filename);
p2f = dsswrite(tc,'ieee13',0,'ieee13',cmds);

[a2 b2] = dssget(p2f);
dssstructcmp(b1,b2);

%% test IEEE bus30 circuit
% run simulation on original file provided in OpenDSS folder
filename = 'IEEETestCases\IEEE 30 Bus\Run_IEEE30.DSS';
[a1 b1] = dssget(filename);

% convert to opendss struct and run it
[tc cmds] = dssparse(filename);
p2f = dsswrite(tc,'ieee30',0,'ieee30',cmds);

[a2 b2] = dssget(p2f);
dssstructcmp(b1,b2);

%% test IEEE bus34 circuit
% run simulation on original file provided in OpenDSS folder
filename = 'IEEETestCases\34Bus\Run_IEEE34Mod1.dss';
[a1 b1] = dssget(filename);

% convert to opendss struct and run it
[tc cmds] = dssparse(filename);
p2f = dsswrite(tc,'ieee30',0,'ieee30',cmds);

[a2 b2] = dssget(p2f);
dssstructcmp(b1,b2);

%% test IEEE bus37 circuit
% run simulation on original file provided in OpenDSS folder
filename = 'IEEETestCases\IEEE 30 Bus\Run_IEEE30.DSS';
[a1 b1] = dssget(filename);

% convert to opendss struct and run it
[tc cmds] = dssparse(filename);
p2f = dsswrite(tc,'ieee30',0,'ieee30',cmds);

[a2 b2] = dssget(p2f);
dssstructcmp(b1,b2);

%% These should be working but you might consider if you really want to run
% them. They take a while.
% Test EPRITestCase ckt24
[c24 cmd24] = dssparse('EPRITestCircuits\ckt24/Run_Ckt24.dss');
p2f = dsswrite(c24,'ckt24',1,'ckt24',cmd24);

dssget('EPRITestCircuits\ckt24/Run_Ckt24.dss');
dssget(p2f);

%% Test EPRITestCase ckt7
[c7 cmd7] = dssparse('EPRITestCircuits\ckt7/RunDSS_Ckt7.dss');
p2f = dsswrite(c7,'ckt7',1,'ckt7',cmd7);

dssget('EPRITestCircuits\ckt7/RunDSS_Ckt7.dss');
dssget(p2f);

%% Test EPRITestCase ckt5
[c5 cmd5] = dssparse('EPRITestCircuits\ckt5/Run_Ckt5.dss');
p2f = dsswrite(c5,'ckt5',1,'ckt5',cmd5);

dssget('EPRITestCircuits\ckt5/Run_Ckt5.dss');
dssget(p2f);

%% TODO: does not work on following circuits yet
%% Test 8500 node feeder
[c85 cmd85] = dssparse('IEEETestCases/8500-Node/Run_8500Node.dss');

%% No need to test following parts
%% TEST: Find unknown linecodes
% Line resistance unit: ohm/1000ft
dat = {d1, d2, d3, d4, d5};
unknownCond = [];
for k = 1:length(dat)
    d = dat{k};
    lc = excel2obj( [curdir '\linecode.xlsx'] );

    % Check for missing linecode
    cond = unique([d.InstSection.PhaseConductorId; d.InstSection.NeutralConductorId]);
    fcond = lc.LineCode.ConductorId;

    id = [];
    for i = 1:length(cond)
        if sum( strcmp(cond{i}, fcond) ) < 1
            id = [id i];
        end
    end
    unknownCond = unique([unknownCond; cond(id)]);
end

%% TEST: Draw Loads' position
close all
drawSections('generate',d);
drawSections(d.InstSection.SectionId,'b');
drawSections(d.Loads.SectionId,'r',1);
hold on;

% node indexing based on names
% build a hash lookup table for node coords
	for i=1:length(d.Node.NodeId)
		nodeXY.(['n' d.Node.NodeId{i}]) = [d.Node.X{i} d.Node.Y{i}];
	end

% section indexing based on names
% build a hash lookup table for sections
for i=1:length(d.InstSection.SectionId)
    id = ['s' d.InstSection.SectionId{i}];
    id(id==' ') = [];
    sectionH.(id) = struct('FromNodeId',d.InstSection.FromNodeId{i},'ToNodeId',d.InstSection.ToNodeId{i});
end
%
for i=1:length(d.Loads.SectionId)
	id = ['s' d.Loads.SectionId{i}];
	id(id==' ') = [];
	section = sectionH.(id);
	coords = vertcat(nodeXY.(['n' section.ToNodeId]));
	plot(coords(1),coords(2),'o');
end