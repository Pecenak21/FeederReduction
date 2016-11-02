%% fix the phases for PVsystems
load data/f520.mat;
pvPro = load('Fallbrook_Scenario2/scenario2_pvProfiles.mat');
pv = load('data/scenario2_pv.mat'); pv = pv.pv; 

%% find 3-phase big loads
l3_id = find([c.load.Phases] == 3);

%% find 3-phase pv (exclude 45 original systems)
p3_id = find([pv.phases] == 3);
p3_id = setdiff(p3_id, [1:45]);

%% randomly assign 3-phase pv to 3-phase loads
% save new bus (DO NOT RUN THIS BLOCK since we don't want to regenerate random buses)
% bus_id = randperm(length(l3_id));
% bus_id = bus_id(1:length(p3_id));
% bus_id = l3_id(bus_id);
% save('data/bus3phase.mat','l3_id','bus_id','p3_id');

%% assigning bus...
load('data/bus3phase.mat');
for i = 1:length(p3_id)
	pv(p3_id(i)).bus1 = c.load(bus_id(i)).bus1;
end

%% compare PV rated and load rated power
pp = pv(p3_id).kVA;
lp = c.load(bus_id).Kw;
co = zeros(length(pp),2);
for i = 1:length(pp)
	co(i,:) = [pp{i} lp{i}];
end

%% save pv system
save('data/scenario2_pv.mat','pv');

%% check location of all 3-phase pv sites
c2 = c;
c2.pvsystem = pv(p3_id);
circuitVisualizer(c2)

%% PV systems that have larger rating than loads
id = find(co(:,1) > co(:,2));
c3 = c;
c3.pvsystem = pv(p3_id(id));
circuitVisualizer(c3);