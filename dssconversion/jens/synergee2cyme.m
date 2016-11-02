function [ ob ] = synergee2cyme(o)
% synergee2cyme( o)
%
% PURPOSE : Converts object from synergee data to object in cyme format.
%           The cyme-formated object can be used as input to cyme2obj fnct.
%
%
% INPUT :   o: object that is created from synergee data by excel2obj fnct.
%
%
% OUTPUT :  a structure that follows cyme format


InstSection=structconv(o.InstSection);
Node=structconv(o.Node);
InstFeeders=structconv(o.InstFeeders);
Loads=structconv(o.Loads);
if isfield(o,'InstLargeCust')
    LargeCust=structconv(o.InstLargeCust);
end
glc=structconv(o.LineCode);
SAI_Control=structconv(o.SAI_Control);
if isfield(o,'InstFuses')
    InstFuses=structconv(o.InstFuses);
end
if isfield(o,'InstCapacitors')
    InstCapacitors=structconv(o.InstCapacitors);
end
if isfield(o,'InstPrimaryTransformers')
    InstTransformers=structconv(o.InstPrimaryTransformers); % Currently not being used as there are no primary transformers specified in SDGE system
end
if isfield(o,'InstRegulators')
    InstRegulators=structconv(o.InstRegulators);
end
if isfield(o,'InstSwitches')
    InstSwitches=structconv(o.InstSwitches);
end
if isfield(o,'InstReclosers')
    InstReclosers=structconv(o.InstReclosers);
end


%%  populate ob.network.switch_setting

switch_setting.SectionID = ...
    InstSwitches.SectionId;

switch_setting.DeviceNumber = ...
    InstSwitches.UniqueDeviceId;

switch_setting.EqID = ...
    InstSwitches.SwitchType;

switch_setting.Automated = ...
    InstSwitches.IsAutomaticSwitch;

switch_setting.SwitchIsOpen = ... % Only used in SynerGEE
    InstSwitches.SwitchIsOpen;

switch_setting.ClosedPhase=num2cellstr(cell2mat(InstSwitches.SwitchIsOpen));

index_SwitchIsOpen=ismember(switch_setting.ClosedPhase,{'1'});

switch_setting.ClosedPhase(index_SwitchIsOpen)={'NONE'};
switch_setting.ClosedPhase(~index_SwitchIsOpen)={'ABC'}; % This assumes that either all phases are open, or all are closed

ob.network.switch_setting=structconv(switch_setting);


%%  populate ob.network.recloser_setting

recloser_setting.SectionID = ...
    InstReclosers.SectionId;

recloser_setting.DeviceNumber = ...
    InstReclosers.UniqueDeviceId;

recloser_setting.NearFromNode = ... % Only used in SynerGEE
    InstReclosers.NearFromNode;

recloser_setting.FastPhaseCount = ...  % Only used in SynerGEE
    InstReclosers.FastPhaseCount;

recloser_setting.SlowPhaseCount = ...  % Only used in SynerGEE
    InstReclosers.SlowPhaseCount;

recloser_setting.SlowGroundTimeAddSec = ... % Only used in SynerGEE
    InstReclosers.SlowGroundTimeAddSec;

recloser_setting.FastGroundTimeAddSec = ... % Only used in SynerGEE
    InstReclosers.FastGroundTimeAddSec;

recloser_setting.SlowPhaseTimeAddSec = ... % Only used in SynerGEE
    InstReclosers.SlowPhaseTimeAddSec;

recloser_setting.FastPhaseTimeAddSec = ... % Only used in SynerGEE
    InstReclosers.FastPhaseTimeAddSec;

recloser_setting.ClosedPhase=num2cellstr(cell2mat(InstReclosers.RecloserIsOpen));

index_RecloserIsOpen=ismember(recloser_setting.ClosedPhase,{'1'});

recloser_setting.ClosedPhase(index_RecloserIsOpen)={'NONE'};
recloser_setting.ClosedPhase(~index_RecloserIsOpen)={'ABC'}; % This assumes that either all phases are open, or all are closed

ob.network.recloser_setting=structconv(recloser_setting);

%%  populate ob.network.section.section

section.SectionID = ...
    InstSection.SectionId;

section.NetworkID = ...
    InstSection.FeederId;
     
section.FromNodeID = ...
    InstSection.FromNodeId;

section.ToNodeID = ...
    InstSection.ToNodeId;

section.Phase = ...
    InstSection.SectionPhases;

section.Phase1 = ...
    SAI_Control.Phase1;

section.Phase2 = ...
    SAI_Control.Phase2;

section.Phase3 = ...
    SAI_Control.Phase3;

ob.network.section.section=structconv(section);

%%  populate ob.network.section.feeder (may not be necessary)

feeder.NetworkID = ...
    InstFeeders.SubstationId;

feeder.HeadNodeID = ...
    InstFeeders.FeederId;

ob.network.section.feeder=structconv(feeder);

%%  populate ob.network.node

node.NodeID = ...
    Node.NodeId;

node.CoordX = ...
    Node.X;

node.CoordY = ...
    Node.Y;

ob.network.node=structconv(node);

%%  populate ob.network.source

source.NodeID = ...
    InstFeeders.FeederId;

source.SourceID = ...
    InstFeeders.FeederId;

source.NetworkID = ...
    InstFeeders.SubstationId;

ob.network.source=structconv(source);


%%  populate ob.network.source_equivalent
% data redundant with ob.equipment.substation

source_equivalent.NodeID = ...
    InstFeeders.FeederId;

source_equivalent.Voltage = ...
    num2cellstr(cell2mat(InstFeeders.NominalKvll));

source_equivalent.PositiveSequenceResistance = ...
    num2cellstr(cell2mat(InstFeeders.PosSequenceResistance));

source_equivalent.PositiveSequenceReactance = ...
    num2cellstr(cell2mat(InstFeeders.PosSequenceReactance));

source_equivalent.ZeroSequenceResistance = ...
    num2cellstr(cell2mat(InstFeeders.ZeroSequenceResistance));

source_equivalent.ZeroSequenceReactance = ...
    num2cellstr(cell2mat(InstFeeders.ZeroSequenceReactance));



% substation.ImpedanceUnit = ...;

ob.network.source_equivalent=structconv(source_equivalent);


%%  populate ob.equipment.substation
% not sure where MVA is specified in synergee

substation.ID = ...
    InstFeeders.FeederId;

substation.KVLL = ...
    num2cellstr(cell2mat(InstFeeders.NominalKvll));

substation.R1 = ...
    num2cellstr(cell2mat(InstFeeders.PosSequenceResistance));

substation.X1 = ...
    num2cellstr(cell2mat(InstFeeders.PosSequenceReactance));

substation.R0 = ...
    num2cellstr(cell2mat(InstFeeders.ZeroSequenceResistance));

substation.X0 = ...
    num2cellstr(cell2mat(InstFeeders.ZeroSequenceReactance));

substation.Conn = ...
    InstFeeders.ConnectionType;

substation.ByPhVoltDegPh1 = ...
    num2cellstr(cell2mat(InstFeeders.ByPhVoltDegPh1));

substation.ByPhVoltDegPh2 = ...
    num2cellstr(cell2mat(InstFeeders.ByPhVoltDegPh2));

ob.equipment.substation=structconv(substation);

%%  populate ob.network.overheadline_setting

overheadline_setting.SectionID = ...
    InstSection.SectionId;

overheadline_setting.LineCableID = ...
    InstSection.PhaseConductorId;
 
overheadline_setting.Length = ...
    InstSection.SectionLength_MUL;

ob.network.overheadline_setting=structconv(overheadline_setting);

%%  populate ob.load.customer_loads

% Phase A
% get indeces of all non-zero entries
Phase1Kw_num=cell2mat(Loads.Phase1Kw);
index_A=find(Phase1Kw_num); % assuming all loaded phases have real power
SectionID1=Loads.SectionId(index_A);
CustomerType1=Loads.SectionId(index_A);
Phase1Kw_num=Phase1Kw_num(index_A);
Phase1Kvar_num=cell2mat(Loads.Phase1Kvar);
Phase1Kvar_num=Phase1Kvar_num(index_A);
Phase1Kva_num=cell2mat(Loads.Phase1Kva);
Phase1Kva_num=Phase1Kva_num(index_A);
Phase1Kwh_num=cell2mat(Loads.Phase1Kwh);
Phase1Kwh_num=Phase1Kwh_num(index_A);
S1_num=sqrt(Phase1Kw_num.^2+Phase1Kvar_num.^2);
pf1_perc_num=Phase1Kw_num./S1_num*100;
Phase1Customers_num=cell2mat(Loads.Phase1Customers);
Phase1Customers_num=Phase1Customers_num(index_A);

Phase1Kw=num2cellstr(Phase1Kw_num);
Phase1Kvar=num2cellstr(Phase1Kvar_num);
Phase1Kva=num2cellstr(Phase1Kva_num);
Phase1Kwh=num2cellstr(Phase1Kwh_num);
pf1_perc=num2cellstr(pf1_perc_num);
S1=num2cellstr(S1_num);
PhaseA(1:length(index_A),1)=deal({'A'});
Phase1Customers=num2cellstr(Phase1Customers_num);


% Phase B
% get indeces of all non-zero entries
Phase2Kw_num=cell2mat(Loads.Phase2Kw);
index_B=find(Phase2Kw_num); % assuming all loaded phases have real power
SectionID2=Loads.SectionId(index_B);
CustomerType2=Loads.SectionId(index_B);
Phase2Kw_num=Phase2Kw_num(index_B);
Phase2Kvar_num=cell2mat(Loads.Phase2Kvar);
Phase2Kvar_num=Phase2Kvar_num(index_B);
Phase2Kva_num=cell2mat(Loads.Phase2Kva);
Phase2Kva_num=Phase2Kva_num(index_B);
Phase2Kwh_num=cell2mat(Loads.Phase2Kwh);
Phase2Kwh_num=Phase2Kwh_num(index_B);
S2_num=sqrt(Phase2Kw_num.^2+Phase2Kvar_num.^2);
pf2_perc_num=Phase2Kw_num./S2_num*100;
Phase2Customers_num=cell2mat(Loads.Phase2Customers);
Phase2Customers_num=Phase2Customers_num(index_B);

Phase2Kw=num2cellstr(Phase2Kw_num);
Phase2Kvar=num2cellstr(Phase2Kvar_num);
Phase2Kva=num2cellstr(Phase2Kva_num);
Phase2Kwh=num2cellstr(Phase2Kwh_num);
pf2_perc=num2cellstr(pf2_perc_num);
S2=num2cellstr(S2_num);
PhaseB(1:length(index_B),1)=deal({'B'});
Phase2Customers=num2cellstr(Phase2Customers_num);


% Phase C
% get indeces of all non-zero entries
Phase3Kw_num=cell2mat(Loads.Phase3Kw);
index_C=find(Phase3Kw_num); % assuming all loaded phases have real power
SectionID3=Loads.SectionId(index_C);
CustomerType3=Loads.SectionId(index_C);
Phase3Kw_num=Phase3Kw_num(index_C);
Phase3Kvar_num=cell2mat(Loads.Phase3Kvar);
Phase3Kvar_num=Phase3Kvar_num(index_C);
Phase3Kva_num=cell2mat(Loads.Phase3Kva);
Phase3Kva_num=Phase3Kva_num(index_C);
Phase3Kwh_num=cell2mat(Loads.Phase3Kwh);
Phase3Kwh_num=Phase3Kwh_num(index_C);
S3_num=sqrt(Phase3Kw_num.^2+Phase3Kvar_num.^2);
pf3_perc_num=Phase3Kw_num./S3_num*100;
Phase3Customers_num=cell2mat(Loads.Phase3Customers);
Phase3Customers_num=Phase3Customers_num(index_C);

Phase3Kw=num2cellstr(Phase3Kw_num);
Phase3Kvar=num2cellstr(Phase3Kvar_num);
Phase3Kva=num2cellstr(Phase3Kva_num);
Phase3Kwh=num2cellstr(Phase3Kwh_num);
pf3_perc=num2cellstr(pf3_perc_num);
S3=num2cellstr(S3_num);
PhaseC(1:length(index_C),1)=deal({'C'});
Phase3Customers=num2cellstr(Phase3Customers_num);

customer_loads.SectionID=[SectionID1;SectionID2;SectionID3];
customer_loads.CustomerType=[CustomerType1;CustomerType2;CustomerType3];
customer_loads.Value1=[Phase1Kw;Phase2Kw;Phase3Kw];
customer_loads.Value2=[pf1_perc;pf2_perc;pf3_perc];
customer_loads.S=[S1;S2;S3];
customer_loads.Kvar=[Phase1Kvar;Phase2Kvar;Phase3Kvar];
customer_loads.ConnectedKVA=[Phase1Kva;Phase2Kva;Phase3Kva];
customer_loads.KWH=[Phase1Kwh;Phase2Kwh;Phase3Kwh];
customer_loads.LoadPhase=[PhaseA;PhaseB;PhaseC];
customer_loads.NumberOfCustomer=[Phase1Customers;Phase2Customers;Phase3Customers];

% Assuming all loads are spot loads, even if they are not. I am making this
% assumptions because OpenDSS cannot handle distributed loads.
LoadType(1:length([index_A;index_B;index_C]),1)=deal({'SPOT'});
customer_loads.LoadType=LoadType;

LocationToModelSpotLoads(1:length([index_A;index_B;index_C]),1)=deal(Loads.LocationToModelSpotLoads(1)); % assuming that the first entry and all other entries are identical, 'F' means load at from node, 'C' means load at to node
customer_loads.LocationToModelSpotLoads=LocationToModelSpotLoads;


ob.load.customer_loads=structconv(customer_loads);

%%  populate ob.load.large_customer_loads
% This object does not exist in CYME. Large customers have three phase
% loads and some of them have generation. The dssconversion_CYME function
% processes this information if the large_customer object exists and
% creates loads and generators.

% Phase A
% get indeces of all non-zero entries
Phase1Kw_num=cell2mat(LargeCust.LoadPhase1Kw);
index_A=find(Phase1Kw_num); % assuming all loaded phases have real power
SectionID1=LargeCust.SectionId(index_A);
CustomerType1=LargeCust.SectionId(index_A);
Phase1Kw_num=Phase1Kw_num(index_A);
Phase1Kvar_num=cell2mat(LargeCust.LoadPhase1Kvar);
Phase1Kvar_num=Phase1Kvar_num(index_A);
Phase1Kva_num=cell2mat(LargeCust.Phase1kva);
Phase1Kva_num=Phase1Kva_num(index_A);
% Phase1Kwh_num=cell2mat(LargeCust.Phase1Kwh);
% Phase1Kwh_num=Phase1Kwh_num(index_A);
S1_num=sqrt(Phase1Kw_num.^2+Phase1Kvar_num.^2);
pf1_perc_num=Phase1Kw_num./S1_num*100;
Phase1Customers_num=cell2mat(LargeCust.Phase1Customers);
Phase1Customers_num=Phase1Customers_num(index_A);

Phase1Kw=num2cellstr(Phase1Kw_num);
Phase1Kvar=num2cellstr(Phase1Kvar_num);
Phase1Kva=num2cellstr(Phase1Kva_num);
% Phase1Kwh=num2cellstr(Phase1Kwh_num);
pf1_perc=num2cellstr(pf1_perc_num);
S1=num2cellstr(S1_num);
clear PhaseA;
PhaseA(1:length(index_A),1)=deal({'A'});
Phase1Customers=num2cellstr(Phase1Customers_num);

% Phase B
% get indeces of all non-zero entries
Phase2Kw_num=cell2mat(LargeCust.LoadPhase2Kw);
index_B=find(Phase2Kw_num); % assuming all loaded phases have real power
SectionID2=LargeCust.SectionId(index_B);
CustomerType2=LargeCust.SectionId(index_B);
Phase2Kw_num=Phase2Kw_num(index_B);
Phase2Kvar_num=cell2mat(LargeCust.LoadPhase2Kvar);
Phase2Kvar_num=Phase2Kvar_num(index_B);
Phase2Kva_num=cell2mat(LargeCust.Phase2kva);
Phase2Kva_num=Phase2Kva_num(index_B);
% Phase2Kwh_num=cell2mat(LargeCust.Phase2Kwh);
% Phase2Kwh_num=Phase2Kwh_num(index_B);
S2_num=sqrt(Phase2Kw_num.^2+Phase2Kvar_num.^2);
pf2_perc_num=Phase2Kw_num./S2_num*100;
Phase2Customers_num=cell2mat(LargeCust.Phase2Customers);
Phase2Customers_num=Phase2Customers_num(index_B);

Phase2Kw=num2cellstr(Phase2Kw_num);
Phase2Kvar=num2cellstr(Phase2Kvar_num);
Phase2Kva=num2cellstr(Phase2Kva_num);
% Phase2Kwh=num2cellstr(Phase2Kwh_num);
pf2_perc=num2cellstr(pf2_perc_num);
S2=num2cellstr(S2_num);
clear PhaseB;
PhaseB(1:length(index_B),1)=deal({'B'});
Phase2Customers=num2cellstr(Phase2Customers_num);

% Phase C
% get indeces of all non-zero entries
Phase3Kw_num=cell2mat(LargeCust.LoadPhase3Kw);
index_C=find(Phase3Kw_num); % assuming all loaded phases have real power
SectionID3=LargeCust.SectionId(index_C);
CustomerType3=LargeCust.SectionId(index_C);
Phase3Kw_num=Phase3Kw_num(index_C);
Phase3Kvar_num=cell2mat(LargeCust.LoadPhase3Kvar);
Phase3Kvar_num=Phase3Kvar_num(index_C);
Phase3Kva_num=cell2mat(LargeCust.Phase3kva);
Phase3Kva_num=Phase3Kva_num(index_C);
% Phase3Kwh_num=cell2mat(LargeCust.Phase3Kwh);
% Phase3Kwh_num=Phase3Kwh_num(index_C);
S3_num=sqrt(Phase3Kw_num.^2+Phase3Kvar_num.^2);
pf3_perc_num=Phase3Kw_num./S3_num*100;
Phase3Customers_num=cell2mat(LargeCust.Phase3Customers);
Phase3Customers_num=Phase3Customers_num(index_C);

Phase3Kw=num2cellstr(Phase3Kw_num);
Phase3Kvar=num2cellstr(Phase3Kvar_num);
Phase3Kva=num2cellstr(Phase3Kva_num);
% Phase3Kwh=num2cellstr(Phase3Kwh_num);
pf3_perc=num2cellstr(pf3_perc_num);
S3=num2cellstr(S3_num);
clear PhaseC;
PhaseC(1:length(index_C),1)=deal({'C'});
Phase3Customers=num2cellstr(Phase3Customers_num);

large_customer_loads.SectionID=[SectionID1;SectionID2;SectionID3];
large_customer_loads.CustomerType=[CustomerType1;CustomerType2;CustomerType3];
large_customer_loads.Value1=[Phase1Kw;Phase2Kw;Phase3Kw];
large_customer_loads.Value2=[pf1_perc;pf2_perc;pf3_perc];
large_customer_loads.S=[S1;S2;S3];
large_customer_loads.Kvar=[Phase1Kvar;Phase2Kvar;Phase3Kvar];
large_customer_loads.ConnectedKVA=[Phase1Kva;Phase2Kva;Phase3Kva];
% large_customer_loads.KWH=[Phase1Kwh;Phase2Kwh;Phase3Kwh];
large_customer_loads.LoadPhase=[PhaseA;PhaseB;PhaseC];
large_customer_loads.NumberOfCustomer=[Phase1Customers;Phase2Customers;Phase3Customers];

% Assuming all loads are spot loads, even if they are not. I am making this
% assumptions because OpenDSS cannot handle distributed loads.
clear LoadType;
LoadType(1:length([index_A;index_B;index_C]),1)=deal({'SPOT'});
large_customer_loads.LoadType=LoadType;

clear LocationToModelSpotLoads;
LocationToModelSpotLoads(1:length([index_A;index_B;index_C]),1)=deal(Loads.LocationToModelSpotLoads(1)); % assuming that the first entry and all other entries are identical, 'F' means load at from node, 'C' means load at to node
large_customer_loads.LocationToModelSpotLoads=LocationToModelSpotLoads;


% This part is different from the 'Loads' struct. The difference is that
% the ZIP information for the 'Loads' struct resides in the 'InstSection'
% struct. The 'InstLargeCust' has ZIP information.
LoadPctConstZ=[LargeCust.LoadPctConstZ;LargeCust.LoadPctConstZ;LargeCust.LoadPctConstZ];
LoadPctConstI=[LargeCust.LoadPctConstI;LargeCust.LoadPctConstI;LargeCust.LoadPctConstI];
LoadPctConstZ_num=cell2mat(LoadPctConstZ);
LoadPctConstI_num=cell2mat(LoadPctConstI);
LoadPctConstP_num=ones(length(LoadPctConstZ_num),1)*100-LoadPctConstZ_num-LoadPctConstI_num;
% LoadPctConstP=num2cellstr(LoadPctConstP_num);
large_customer_loads.ConstantImpedance=num2cellstr(LoadPctConstZ_num);
large_customer_loads.ConstantCurrent=num2cellstr(LoadPctConstI_num);
large_customer_loads.ConstantPower=num2cellstr(LoadPctConstP_num);

large_customer_loads.GenTotalKw=num2cellstr(cell2mat([LargeCust.GenTotalKw;LargeCust.GenTotalKw;LargeCust.GenTotalKw]));
large_customer_loads.GenTotalKvar=num2cellstr(cell2mat([LargeCust.GenTotalKvar;LargeCust.GenTotalKvar;LargeCust.GenTotalKvar]));
large_customer_loads.GenStatus=[LargeCust.GenStatus;LargeCust.GenStatus;LargeCust.GenStatus];
large_customer_loads.GenPct=num2cellstr(cell2mat([LargeCust.GenPct;LargeCust.GenPct;LargeCust.GenPct]));
large_customer_loads.GenCustClass=[LargeCust.GenCustClass;LargeCust.GenCustClass;LargeCust.GenCustClass];
large_customer_loads.UniqueDeviceId=[LargeCust.UniqueDeviceId;LargeCust.UniqueDeviceId;LargeCust.UniqueDeviceId];

ob.load.large_customer_loads=structconv(large_customer_loads);

%%  populate ob.load.customer_class
% maped to customer_loads via CustomerType, in Synergee, CustomerType info
% is in InstSection and associated with a section, i.e., each section is 
% CustomerType resulting in lot's of redundancy, could be reduced to one
% (or a few) CustomerTypes to increase speed

customer_class.ID = ...
    InstSection.SectionId;

customer_class.ConstantPower = ...
    InstSection.FeederId;

PercentSpotLoadConstCurrent=cell2mat(InstSection.PercentSpotLoadConstCurrent);
PercentSpotLoadConstImpedance=cell2mat(InstSection.PercentSpotLoadConstImpedance);
PercentSpotLoadConstPower=100-PercentSpotLoadConstCurrent-PercentSpotLoadConstImpedance;

PercentDistLoadConstCurrent=cell2mat(InstSection.PercentDistLoadConstCurrent);
PercentDistLoadConstImpedance=cell2mat(InstSection.PercentDistLoadConstImpedance);
PercentDistLoadConstPower=100-PercentDistLoadConstCurrent-PercentDistLoadConstImpedance;

% using the distributed load characteristics, because, in the SDGE system,
% almost all loads are distributed loads (even though we model them as spot
% loads due to OpenDSS constraints
customer_class.ConstantCurrent = ...
    num2cellstr(PercentDistLoadConstCurrent);

customer_class.ConstantImpedance = ...
    num2cellstr(PercentDistLoadConstImpedance);

customer_class.ConstantPower = ...
    num2cellstr(PercentDistLoadConstPower);


ob.load.customer_class=structconv(customer_class);

%%  populate ob.LineCode

LineCode.ID = ...
    glc.ConductorId;

LineCode.R1 = ...
    num2cellstr(cell2mat(glc.R1));

LineCode.X1 = ...
    num2cellstr(cell2mat(glc.X1));

LineCode.B1 = ...
    num2cellstr(cell2mat(glc.B1));

LineCode.Z1 = ...
    num2cellstr(cell2mat(glc.Z1));

LineCode.R0 = ...
    num2cellstr(cell2mat(glc.R0));

LineCode.X0 = ...
    num2cellstr(cell2mat(glc.X0));

LineCode.B0 = ...
    num2cellstr(cell2mat(glc.B0));

LineCode.Z0 = ...
    num2cellstr(cell2mat(glc.Z0));

LineCode.Rs = ...
    num2cellstr(cell2mat(glc.Rs));

LineCode.Xs = ...
    num2cellstr(cell2mat(glc.Xs));

LineCode.Zs = ...
    num2cellstr(cell2mat(glc.Zs));


ob.LineCode=structconv(LineCode);

%%  populate ob.Fuses

Fuses.SectionID = ...
    InstFuses.SectionId;

Fuses.DeviceNumber = ...
    InstFuses.UniqueDeviceId;

Fuses.Amps = ...
    InstFuses.AmpRating;

index=ismember(num2cellstr(cell2mat(InstFuses.FuseIsOpen)),'0');
InstFuses.FuseIsOpen(index)={'ABC'}; % Ignoring the possibility that a fuse might be on a single-phase section
InstFuses.FuseIsOpen(~index)={''}; % Fuse is open (Fuse should be closed during normal operation)
Fuses.ClosedPhase=InstFuses.FuseIsOpen;

clear index;
index=ismember(num2cellstr(cell2mat(InstFuses.NearFromNode)),'1');
InstFuses.NearFromNode(index)={'L'}; % assuming 'L' means 'NearFromNode' (building the system from left to right), not sure if this assumption is accurate
InstFuses.NearFromNode(~index)={'R'}; % Not sure if 'R' is used in CYME to designate 'not near from node'
Fuses.Location=InstFuses.NearFromNode;
ob.Fuses=Fuses;


%%  populate ob.Capacitors

Capacitors.SectionID = ...
    InstCapacitors.SectionId;

Capacitors.DeviceNumber = ...
    InstCapacitors.UniqueDeviceId;

Capacitors.Phases_xyz = ...
    InstCapacitors.ConnectedPhases; 

Capacitors.Phases = ...
    num2cell(cellfun(@length,InstCapacitors.ConnectedPhases));

Capacitors.Connection = ...
    InstCapacitors.ConnectionType;

Capacitors.ConnectionStatus = ...
    InstCapacitors.CapacitorIsOn;


Capacitors.Module1On = ...
    InstCapacitors.Module1On;
Capacitors.Module1Activated = ...
    InstCapacitors.Module1Activated;
Capacitors.Module1KvarPerPhase = ...
    InstCapacitors.Module1KvarPerPhase;
Capacitors.Module1CapSwitchCloseValue = ...
    InstCapacitors.Module1CapSwitchCloseValue;
Capacitors.Module1CapSwitchTripValue = ...
    InstCapacitors.Module1CapSwitchTripValue;

Capacitors.Module2On = ...
    InstCapacitors.Module2On;
Capacitors.Module2Activated = ...
    InstCapacitors.Module2Activated;
Capacitors.Module2KvarPerPhase = ...
    InstCapacitors.Module2KvarPerPhase;
Capacitors.Module2CapSwitchCloseValue = ...
    InstCapacitors.Module2CapSwitchCloseValue;
Capacitors.Module2CapSwitchTripValue = ...
    InstCapacitors.Module2CapSwitchTripValue;

Capacitors.Module3On = ...
    InstCapacitors.Module3On;
Capacitors.Module3Activated = ...
    InstCapacitors.Module3Activated;
Capacitors.Module3KvarPerPhase = ...
    InstCapacitors.Module3KvarPerPhase;
Capacitors.Module3CapSwitchCloseValue = ...
    InstCapacitors.Module3CapSwitchCloseValue;
Capacitors.Module3CapSwitchTripValue = ...
    InstCapacitors.Module3CapSwitchTripValue;

Capacitors.KVLN = ...
    num2cellstr(cell2mat(InstCapacitors.RatedKv));


% SynerGEE specific parameters (not sure what the CYME equivalents are)
Capacitors.PrimaryControlMode = ...
    InstCapacitors.PrimaryControlMode;

Capacitors.MeteringPhase = ...
    InstCapacitors.MeteringPhase;

Capacitors.Module1CapSwitchTripValue = ...
    InstCapacitors.Module1CapSwitchTripValue;
Capacitors.Module2CapSwitchTripValue = ...
    InstCapacitors.Module2CapSwitchTripValue;
Capacitors.Module3CapSwitchTripValue = ...
    InstCapacitors.Module3CapSwitchTripValue;

Capacitors.CapVoltageOverrideActive = ...
    InstCapacitors.CapVoltageOverrideActive;

Capacitors.CapVoltageOverrideSetting = ...
    InstCapacitors.CapVoltageOverrideSetting;

Capacitors.CapVoltageOverrideBandwidth = ...
    InstCapacitors.CapVoltageOverrideBandwidth;

Capacitors.CapacitorPTRatio = ...
    InstCapacitors.CapacitorPTRatio;

Capacitors.CapacitorCTRating = ...
    InstCapacitors.CapacitorCTRating;

fixedKvar=sum([cell2mat(InstCapacitors.FixedKvarPhase1) cell2mat(InstCapacitors.FixedKvarPhase2) cell2mat(InstCapacitors.FixedKvarPhase3)],2);

for i=1:length(InstCapacitors.UniqueDeviceId)
    kvar = Capacitors.Phases{i}*[InstCapacitors.Module1KvarPerPhase{i} InstCapacitors.Module2KvarPerPhase{i} InstCapacitors.Module3KvarPerPhase{i}];
    ind_z = find(kvar==0,1);
    kvar(ind_z:end) = [];

    if fixedKvar(i) > 0
        if(~isempty(kvar)), warning('dssconversion:capacitor:multiplekvar','Capacitor kvar was specified both as fixed and as switchable, which is invalid'); end
        Capacitors.ThreePhaseKVAR{i,1} = num2str(fixedKvar(i));
    elseif ~isempty(kvar)
        if(length(kvar)>1), Capacitors.Numsteps(i)=length(kvar); end
        Capacitors.ThreePhaseKVAR{i,1} = num2str(kvar);
    else
        error('DSSConversion:Capcitor:kvar','Invalid kvar input');
    end
end



% ob.Capacitors=Capacitors;
ob.network.shunt_capacitor_setting=structconv(Capacitors);


%%  populate ob.Regulators

% simply copy content of InstRegulators
% not really a conversion, but good enough for now as we have not dealt
% with CYME regulators, yet
Regulators=InstRegulators;

% Regulators.SectionID = ...
%     InstRegulators.SectionId;
% 
% Regulators.Name = ...
%     InstRegulators.UniqueDeviceId;
% 
% Regulators.Phases_xyz = ...
%     InstRegulators.ConnectedPhases;
% 
% Regulators.RegulatorType = ...
%     InstRegulators.RegulatorType;
% 
% Regulators.RegulatorIsOn = ...
%     InstRegulators.RegulatorIsOn;
% 
% Regulators.RegulatorIsOn = ...
%     InstRegulators.RegulatorIsOn;
% 
% Regulators.NearFromNode = ...
%     InstRegulators.NearFromNode;
% 
% Regulators.ForwardVoltageSettingPhase1 = ...
%     InstRegulators.ForwardVoltageSettingPhase1;
% 
% Regulators.ForwardVoltageSettingPhase2 = ...
%     InstRegulators.ForwardVoltageSettingPhase2;
% 
% Regulators.ForwardVoltageSettingPhase3 = ...
%     InstRegulators.ForwardVoltageSettingPhase3;

ob.Regulators=Regulators;

end