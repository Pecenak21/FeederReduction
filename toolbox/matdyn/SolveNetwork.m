function [U] = SolveNetwork(Xgen, Pgen, Ly, Uy, Py, gbus, gentype)

% [U] = SolveNetwork(Xgen, Pgen, Ly ,Uy ,Py ,gbus, gentype)
% 
% Solve the network
% 
% INPUTS
% Xgen = state variables of generators
% Pgen = parameters of generators
% Ly ,Uy ,Py = factorised augmented bus admittance matrix
% gbus = generator buses
% gentype = generator models
% 
% OUTPUTS
% U = bus voltages
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Init
ngen = length(gbus);
Igen = zeros(ngen,1);

s=length(Py);

Ig = zeros(s,1);
d = [1:length(gentype)]';

%% Define generator types
type1 = d(gentype==1);
type2 = d(gentype==2);

%% Generator type 1: classical model
delta = Xgen(type1,1);
Eq_tr = Xgen(type1,3);

xd_tr = Pgen(type1,7);

% Calculate generator currents
Igen(type1) = (Eq_tr.*exp(j.*delta))./(j.*xd_tr);

%% Generator type 2: 4th order model
delta = Xgen(type2,1);
Eq_tr = Xgen(type2,3);
Ed_tr = Xgen(type2,4);

xd_tr = Pgen(type2,8);

% Calculate generator currents
Igen(type2) = (Eq_tr + j.*Ed_tr).*exp(j.*delta)./(j.*xd_tr);% Padiyar, p.417.

%% Calculations
% Generator currents
Ig(gbus) = Igen;

% Calculate network voltages: U = Y/Ig
U = Uy\(Ly\Ig(Py));

return;