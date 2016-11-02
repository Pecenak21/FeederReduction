clc
clear all

load('C:\Work\Projects\2012\1787-UCSD PV\Simulation\System\SDGE\tmp/f520.mat');
c.XFR_Load = c.transformer;
% insert single phase trafos servicing single phase loads
loads_1Phase = c.load([c.load.Phases]==1);
[a0 b0] = ismember({loads_1Phase.Name}, {c.load.Name});
for j=1:1:length(b0)
    nl = dsstransformer('Name',strcat('t_', c.load(b0(j)).Name));
    nl.Phases = c.load(b0(j)).Phases;
    nl.kVs = [12.47/sqrt(3) 0.240];%For 2-or 3-phase, enter phase-phase kV rating.  Otherwise, kV rating of the actual winding. 240V for the sec corresponds to single phase 3 wire connection.
    KVAs = 2*sqrt(c.load(b0(j)).Kw^2 + c.load(b0(j)).Kvar^2);%rate the trafo at twice the load's actual kVA (not default kVA)
    nl.kVAs = [KVAs KVAs];
    nl.imag = 0.5;%default value from DG system
    nl.Noloadloss = 0.25;%default value from DG system
    LoadBus = strcat('X_',c.load(b0(j)).bus1);%change trafo sec bus and load bus names
    nl.Buses = {c.load(b0(j)).bus1 LoadBus};
    c.load(b0(j)).bus1 = LoadBus;
    c.load(b0(j)).Kv = 0.240;% set single phase loads LN voltage levels to 240V
    c.XFR_Load(end+1)=nl;
end
% insert three phase trafos servicing three phase loads
loads_3Phase = c.load([c.load.Phases]==3);
[a0 b0] = ismember({loads_3Phase.Name}, {c.load.Name});
for j=1:1:length(b0)
nl = dsstransformer('Name',strcat('t_', c.load(b0(j)).Name));
    nl.Phases = c.load(b0(j)).Phases;
    nl.kVs = [12.47 0.208];%For 2-or 3-phase, enter phase-phase kV rating.  Otherwise, kV rating of the actual winding
    KVAs = 2*sqrt(c.load(b0(j)).Kw^2 + c.load(b0(j)).Kvar^2);%rate the trafo at twice the load's actual kVA (not default kVA)
    nl.kVAs = [KVAs KVAs];
    nl.imag = 0.5;%default value from DG system
    nl.Noloadloss = 0.25;%default value from DG system
    LoadBus = strcat('X_',c.load(b0(j)).bus1);%change trafo sec bus and load bus names
    nl.Buses = {c.load(b0(j)).bus1 LoadBus};
    c.load(b0(j)).bus1 = LoadBus;
    c.load(b0(j)).Kv = 0.208;% set three phase loads LL voltage levels to 208V
    c.XFR_Load(end+1)=nl;
end
bb=1:1:length(c.transformer.Name(:));
c.XFR_Load(bb) = [];
dsswrite(c,'test','1','C:\Work\Projects\2012\1787-UCSD PV\Simulation\System\SDGE\SIMULATION_RUNS\Test_5\WithPV\Test');