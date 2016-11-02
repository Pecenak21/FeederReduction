%% This code is to quasi-state simulation (time series simulation) for a day with deaggregated PV profiles 
clear 
%% Load the circuit with all the modifications in place to match power flow results given from SDG&E
% changes
c = load('tmp/f520_standard.mat','c','glc'); glc = c.glc; c = c.c;

%% run quick power validation check on the circuit
validatepower(3,'kw',0,'',[1 0 0 0 0],c,glc,0);
validatepower(3,'kvar',0,'',[1 0 0 0 0],c,glc,0);