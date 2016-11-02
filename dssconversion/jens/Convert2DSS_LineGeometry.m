function lg = Convert2DSS_LineGeometry(s,glg,name)
% lg = Convert2DSS_LineGeometry(s)
%
% PURPOSE : Converts lines from object to OpenDSS
%
%
% INPUT :   s: struct that has all entries (e.g., o.R1)
%
% OUTPUT :  OpenDSS Lines object
%           see page 106 in OpenDSS manual V7.4.3(March 2012)


lg = dsslinegeometry;
n_phase=0;
n_wire=0;
lg.Name = name;
if ~strcmpi(s.CondID_A,'none')
    lg.Cond=1;
    lg.Wire=s.CondID_A;
    lg.X=glg.PosOfCond1_X;
    lg.H=glg.PosOfCond1_Y;
    lg.Units='m';
    n_phase=n_phase+1;
    n_wire=n_wire+1;
end
if ~strcmpi(s.CondID_B,'none')
    lg.Cond=2;
    lg.Wire=s.CondID_B;
    lg.X=glg.PosOfCond2_X;
    lg.H=glg.PosOfCond2_Y;
    lg.Units='m';
    n_phase=n_phase+1;
    n_wire=n_wire+1;
end
if ~strcmpi(s.CondID_C,'none')
    lg.Cond=3;
    lg.Wire=s.CondID_C;
    lg.X=glg.PosOfCond3_X;
    lg.H=glg.PosOfCond3_Y;
    lg.Units='m';
    n_phase=n_phase+1;
    n_wire=n_wire+1;
end
if ~strcmpi(s.CondID_N,'none')
    lg.Cond=4;
    lg.Wire=s.CondID_N;
    lg.X=glg.PosOfNeutralCond_X;
    lg.H=glg.PosOfNeutralCond_Y;
    lg.Units='m';
    n_wire=n_wire+1;
end
lg.Nconds=n_wire;
lg.Nphases=n_phase;


end



% obj.defaults = struct('Name','', ...
% 		'Nconds',3,... %Number of conductors in this geometry. Default is 3. Triggers memory allocations. Define first!
% 		'Nphases',3,... % Number of phases. Default =3; All other conductors are considered neutrals and might be reduced out.
% 		'Cond',{1},... % Set this to number of the conductor you wish to define. Default is 1.
% 		'Wire',{''},...    % Code from WireData. MUST BE PREVIOUSLY DEFINED. no default.
% 		'X',{[]},... % x coordinate
% 		'H',{[]},... % height of conductor
% 		'Units',{'ft'},... %Units for x and h: {mi|kft|km|m|Ft|in|cm } Initial default is "ft", but defaults to last unit defined
% 		'Normamps',[],... %Normal ampacity, amperes for the line. Defaults to first conductor if not specified.
% 		'Emergamps',[],... %Emergency ampacity, amperes. Defaults to first conductor if not specified.
% 		'Reduce','No',... %{ Yes | No} Default = no. Reduce to Nphases (Kron Reduction). Reduce out neutrals.
% 		'Like','');

% This class of data is used to define the positions of the conductors.
% Nconds= Number of conductors in this geometry. Default is 3. Triggers memory allocations.
% Define first!
% Nphases= Number of phases. Default =3; All other conductors are considered neutrals and
% might be reduced out.
% Cond= Set this to number of the conductor you wish to define. Default is 1.
% Wire= Code from WireData. MUST BE PREVIOUSLY DEFINED. no default.
% X= x coordinate.
% H= Height of conductor.
% Units= Units for x and h: {mi|kft|km|m|Ft|in|cm } Initial default is "ft", but defaults to last
% unit defined
% Normamps= Normal ampacity, amperes for the line. Defaults to first conductor if not specified.
% Emergamps= Emergency ampacity, amperes. Defaults to first conductor if not specified.
% Reduce= { Yes | No} Default = no. Reduce to Nphases (Kron Reduction). Reduce out neutrals.
% Like= Make like another object, e.g.:
% New Capacitor.C2 like=c1 ...
