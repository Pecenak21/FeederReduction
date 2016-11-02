function [ ob ] = cyme2obj( filename,flgChopFile,flgTest,strSystem)

cellSystem_SCE=[{'SCE_Centaur'};{'SCE_Durox'}];
cellSystem_HO=[{'HydroOttawa'};{'default'}];

% Convert cyme files with multiple sheets to object oriented struct
if nargin<3
    flgTest=0;
    flgInput='default';
end
% Chop three main files into smaller files
if(flgChopFile==1)
	GetCYMData_load(filename{1});
    GetCYMData_network(filename{2});
    GetCYMData_equip(filename{3});
end

%% Read in data that resides in smaller files
% - first read data that is in common format
% - next read data that requires individual 'get' functions because of the
%   format of the input files differs
% - Ben is working on 'get' functions that can deal with different formats,
%   which will simplify the code below and should make the conversion of 
%   future CYME systems much easier

% Load data
o.customerloads = GetCYMData_customerloads(filename{1}); % all the load data seems to reside here, the CYMload_LOADS.txt file does not appear to contain any new info
o.customerclass = GetCYMData_customerclass(filename{1});
o.Loads=GetCYMData_loads(filename{1});

% Network data
o.Capacitors=GetCYMData_shuntcapacitorsetting(filename{2}); 
o.recloser = GetCYMData_recloser(filename{2});
o.Node = GetCYMData_node(filename{2});
o.source = GetCYMData_source(filename{2});
o.Fuses = GetCYMData_fusesetting(filename{2});
o.transformersetting = GetCYMData_transformersetting(filename{2});

% Equipment data
o.transformer = GetCYMData_transformer(filename{3});
o.cable = GetCYMData_cable(filename{3});
    

% the sourceequivalent section in the CYME network file
% does not seem to have the correct source information,
% use the info in the substation section (equipment file)instead
o.substation = GetCYMData_substation(filename{3}); 

% this part will become superfluous once we have 'get' functions that can
% deal with differently formatted CYME files
if ismember(strSystem,cellSystem_SCE)
    [o.Section, o.feeder] = GetCYMData_section_SCE(filename{2},flgTest);
    o.sourceequivalent = GetCYMData_sourceequivalent_SCE(filename{2});
    [o.lineconfiguration o.lineconfiguration2] = GetCYMData_LINECONFIGURATION(filename{2});
    o.Switches = GetCYMData_switchsetting_SCE(filename{2});
    o.breaker = GetCYMData_breakersetting_SCE(filename{2});
    o.line = GetCYMData_line_SCE(filename{3});
elseif ismember(strSystem,cellSystem_HO)
    o.Switches = GetCYMData_switchsetting(filename{2});
    o.breaker = GetCYMData_breaker(filename{2});
    [o.Section, o.feeder] = GetCYMData_section(filename{2},flgTest);
    o.undergroundlinesetting = GetCYMData_undergroundlinesetting(filename{2});
    o.sourceequivalent = GetCYMData_sourceequivalent(filename{2});
    o.overheadlinesetting = GetCYMData_overheadlinesetting(filename{2});
    o.line = GetCYMData_line(filename{3});
else
    error('unknown system')
end


%% reduce data
% only include objects and buses that are associated with a section id

% reduce nodes
[flg,index_FromNode]=ismember(o.Section.FromNode,o.Node.NodeID);
[flg,index_ToNode]=ismember(o.Section.ToNode,o.Node.NodeID);
index=[index_ToNode; index_FromNode];
o.Node=structreduce(o.Node,index);

% reduce strucsourceequivalent
[flg,index]=ismember(o.Node.NodeID,o.sourceequivalent.NodeID);
o.sourceequivalent=structreduce(o.sourceequivalent,index);

% reduce strucsource
[flg,index]=ismember(o.Node.NodeID,o.source.NodeID);
o.source=structreduce(o.source,index);

% reduce strucswitchsetting
[flg,index]=ismember(o.Section.SectionID,o.Switches.SectionID);
o.Switches=structreduce(o.Switches,index);

% reduce strucbreaker
[flg,index]=ismember(o.Section.SectionID,o.breaker.SectionID);
o.breaker=structreduce(o.breaker,index);

% reduce strucrecloser
[flg,index]=ismember(o.Section.SectionID,o.recloser.SectionID);
o.recloser=structreduce(o.recloser,index);

% reduce strucfusesetting
[flg,index]=ismember(o.Section.SectionID,o.Fuses.SectionID);
o.Fuses=structreduce(o.Fuses,index);

% reduce structransformersetting
[flg,index]=ismember(o.Section.SectionID,o.transformersetting.SectionID);
o.transformersetting=structreduce(o.transformersetting,index);

% reduce struccustomerloads
[flg,index]=ismember(o.Section.SectionID,o.customerloads.SectionID);
o.customerloads=structreduce(o.customerloads,index);

if ismember(strSystem,cellSystem_SCE)
    % reduce struclineconfiguration
    [flg,index]=ismember(o.Section.SectionID,o.lineconfiguration.SectionID);
    o.lineconfiguration=structreduce(o.lineconfiguration,index);

    % reduce struclineconfiguration2
    [flg,index]=ismember(o.Section.SectionID,o.lineconfiguration2.SectionID);
    o.lineconfiguration2=structreduce(o.lineconfiguration2,index);
elseif ismember(strSystem,cellSystem_HO)
    % reduce strucundergroundlinesetting
    [flg,index]=ismember(o.Section.SectionID,o.undergroundlinesetting.SectionID);
    o.undergroundlinesetting=structreduce(o.undergroundlinesetting,index);

    % reduce strucoverheadlinesetting
    [flg,index]=ismember(o.Section.SectionID,o.overheadlinesetting.SectionID);
    o.overheadlinesetting=structreduce(o.overheadlinesetting,index);
else
    error('unknown system')
end


%% Process Data

% Add information from .transformer to .transformersetting

TSS=o.transformersetting;
SID=TSS.SectionID;
TS=o.transformer; % library of transformers
for i=1:length(SID)
    i_EqID=GetIndex_StrInCell(TS.EqID,TSS.EqID(i));
    o.Transformers.SectionID(i)=TSS.SectionID(i);
    o.Transformers.EqID(i)=TSS.EqID(i);
    o.Transformers.Conn(i)=TSS.Conn(i); %from network data file
    o.Transformers.Conn2(i)=TS.Conn(i_EqID); %from equipment data file
    o.Transformers.kVA(i)=TS.KVA(i_EqID);
    o.Transformers.KVLLprim(i)=TS.KVLLprim(i_EqID);
    o.Transformers.KVLLsec(i)=TS.KVLLsec(i_EqID);
    o.Transformers.Z1(i)=TS.Z1(i_EqID);
    o.Transformers.Z0(i)=TS.Z0(i_EqID);
    o.Transformers.XR(i)=TS.XR(i_EqID);
    o.Transformers.XR0(i)=TS.XR0(i_EqID);
    o.Transformers.Phases(i)={'ABC'}; % assuming transformers are always three phase, not sure where CYME specifies the phasing
    o.Transformers.HighSideNearFromNode(i)=1; % assuming that high side connects to the from node
    connection = regexp(TSS.Conn(i), '_', 'split');
	o.Transformers.HighSideConnectionCode = connection{1}(1);
	o.Transformers.LowSideConnectionCode = connection{1}(2);
    o.Transformers.UniqueDeviceID(i) = {['t_' TSS.SectionID{i}]};
end

% Add information from .customerclass and .customerloads to .Loads
% modify load structure to match format required in dssconversion function
SID=o.Loads.SectionID;
CC=o.customerclass;
CL=o.customerloads;

if ~isempty(o.customerloads.SectionID)
    for i=1:length(SID)
        i_CL=GetIndex_StrInCell(CL.SectionID,SID{i});
        Phase1Kw=0;
        Phase1Value2=0;
        n1Value2=0; % used to calculate average, assuming 'Value2' is power factor
        Phase1Kwh=0;
        Phase1ConnectedKVA=0;
        Phase1Customers=0;

        Phase2Kw=0;
        Phase2Value2=0;
        n2Value2=0;
        Phase2Kwh=0;
        Phase2ConnectedKVA=0;
        Phase2Customers=0;

        Phase3Kw=0;
        Phase3Value2=0;
        n3Value2=0;
        Phase3Kwh=0;
        Phase3ConnectedKVA=0;
        Phase3Customers=0;

        for k=1:length(i_CL) % lump some customer loads together
            if strcmpi(CL.LoadPhase(i_CL(k)),'a')
                Phase1Kw=Phase1Kw+str2double(CL.Value1{i_CL(k)});
                Phase1Value2=Phase1Value2+str2double(CL.Value2{i_CL(k)});
                Phase1Kwh=Phase1Kwh+str2double(CL.KWH{i_CL(k)});
                Phase1ConnectedKVA=Phase1ConnectedKVA+str2double(CL.ConnectedKVA{i_CL(k)});
                Phase1Customers=Phase1Customers+str2double(CL.NumberOfCustomer{i_CL(k)});
                n1Value2=n1Value2+1;
            elseif strcmpi(CL.LoadPhase(i_CL(k)),'b')
                Phase2Kw=Phase2Kw+str2double(CL.Value1{i_CL(k)});
                Phase2Value2=Phase2Value2+str2double(CL.Value2{i_CL(k)});
                Phase2Kwh=Phase2Kwh+str2double(CL.KWH{i_CL(k)});
                Phase2ConnectedKVA=Phase2ConnectedKVA+str2double(CL.ConnectedKVA{i_CL(k)});
                Phase2Customers=Phase2Customers+str2double(CL.NumberOfCustomer{i_CL(k)});
                n2Value2=n2Value2+1;
            elseif strcmpi(CL.LoadPhase(i_CL(k)),'c')
                Phase3Kw=Phase3Kw+str2double(CL.Value1{i_CL(k)});
                Phase3Value2=Phase3Value2+str2double(CL.Value2{i_CL(k)});
                Phase3Kwh=Phase3Kwh+str2double(CL.KWH{i_CL(k)});
                Phase3ConnectedKVA=Phase3ConnectedKVA+str2double(CL.ConnectedKVA{i_CL(k)});
                Phase3Customers=Phase3Customers+str2double(CL.NumberOfCustomer{i_CL(k)});
                n3Value2=n3Value2+1;
            elseif strcmpi(CL.LoadPhase(i_CL(k)),'abc') % not sure if I have to divide the load on each phase by 3
                Phase1Kw=Phase1Kw+str2double(CL.Value1{i_CL(k)});
                Phase1Value2=Phase1Value2+str2double(CL.Value2{i_CL(k)});
                Phase1Kwh=Phase1Kwh+str2double(CL.KWH{i_CL(k)});
                Phase1ConnectedKVA=Phase1ConnectedKVA+str2double(CL.ConnectedKVA{i_CL(k)});
                Phase1Customers=Phase1Customers+str2double(CL.NumberOfCustomer{i_CL(k)});
                n1Value2=n1Value2+1;
                Phase2Kw=Phase2Kw+str2double(CL.Value1{i_CL(k)});
                Phase2Value2=Phase2Value2+str2double(CL.Value2{i_CL(k)});
                Phase2Kwh=Phase2Kwh+str2double(CL.KWH{i_CL(k)});
                Phase2ConnectedKVA=Phase2ConnectedKVA+str2double(CL.ConnectedKVA{i_CL(k)});
                Phase2Customers=Phase2Customers+str2double(CL.NumberOfCustomer{i_CL(k)});
                n2Value2=n2Value2+1;
                Phase3Kw=Phase3Kw+str2double(CL.Value1{i_CL(k)});
                Phase3Value2=Phase3Value2+str2double(CL.Value2{i_CL(k)});
                Phase3Kwh=Phase3Kwh+str2double(CL.KWH{i_CL(k)});
                Phase3ConnectedKVA=Phase3ConnectedKVA+str2double(CL.ConnectedKVA{i_CL(k)});
                Phase3Customers=Phase3Customers+str2double(CL.NumberOfCustomer{i_CL(k)});
                n3Value2=n3Value2+1;
            else
                error('Phase assignment is ' & CC(i_CL(k)) & '. Only a, b, or c allowed (case insensitive).')
            end
        end

            % putting the lumped parameters into a structure
            o.Loads.Phase1Kw{i,1}=Phase1Kw;
            pf1=Phase1Value2/n1Value2/100; % pf=P/S; S=P^2+Q^2; Q=sqrt(P/pf-P^2)
            o.Loads.Phase1Kvar{i,1}=sqrt(Phase1Kw/pf1-Phase1Kw^2);
            o.Loads.Phase1Kwh{i,1}=Phase1Kwh;
            o.Loads.Phase1ConnectedKVA{i,1}=Phase1ConnectedKVA;
            o.Loads.Phase1Customers{i,1}=Phase1Customers;
            o.Loads.Phase2Kw{i,1}=Phase2Kw;
            pf2=Phase2Value2/n2Value2/100;
            o.Loads.Phase2Kvar{i,1}=sqrt(Phase2Kw/pf2-Phase2Kw^2);
            o.Loads.Phase2Kwh{i,1}=Phase2Kwh;
            o.Loads.Phase2ConnectedKVA{i,1}=Phase2ConnectedKVA;
            o.Loads.Phase2Customers{i,1}=Phase2Customers;
            o.Loads.Phase3Kw{i,1}=Phase3Kw;
            pf3=Phase3Value2/n2Value2/100;
            o.Loads.Phase3Kvar{i,1}=sqrt(Phase3Kw/pf3-Phase3Kw^2);
            o.Loads.Phase3Kwh{i,1}=Phase3Kwh;
            o.Loads.Phase3ConnectedKVA{i,1}=Phase3ConnectedKVA;
            o.Loads.Phase3Customers{i,1}=Phase3Customers;

            % assuming that these values do not change for different SectionIDs
            % (or if they change, it does not matter)
            o.Loads.DeviceNumber{i,1}=CL.DeviceNumber(i_CL(1));
            o.Loads.CustomerNumber{i,1}=CL.CustomerNumber(i_CL(1));
            o.Loads.CustomerType{i,1}=CL.CustomerType(i_CL(1));
            o.Loads.ValueType{i,1}=CL.ValueType(i_CL(1));
            o.Loads.ConnectedKVA{i,1}=CL.ConnectedKVA(i_CL(1));
            o.Loads.NumberOfCustomer{i,1}=CL.NumberOfCustomer(i_CL(1));

            i_CC=GetIndex_StrInCell(CC.CustomerType,CL.CustomerType(k)); % i_CC should be unique
            o.Loads.ConstantPower{i,1}=str2double(CC.ConstantPower(i_CC));
            o.Loads.ConstantCurrent{i,1}=str2double(CC.ConstantCurrent(i_CC));
            o.Loads.ConstantImpedance{i,1}=str2double(CC.ConstantImpedance(i_CC));
            o.Loads.UtilizationFactor{i,1}=CC.UtilizationFactor(i_CC);
            o.Loads.PowerFactor{i,1}=CC.PowerFactor(i_CC);
            o.Loads.LoadFactor{i,1}=CC.LoadFactor(i_CC);
            o.Loads.LoadFlowVoltagePercentOfNominal{i,1}=CC.LoadFlowVoltagePercentOfNominal(i_CC);

            % 'F' = load is modeled at the 'from' node, otherwise load is
            % modeled at the 'to' noed
            o.Loads.LocationToModelSpotLoads{i,1}='T';
            o.Loads.LoadType{i,1}=CL.LoadType(i_CL(1));
            if strcmpi(CL.LoadType(i_CL(1)),'spot')
                 o.Loads.IsSpotLoad{i,1}=1;
            else
                o.Loads.IsSpotLoad{i,1}=0;
            end

    end
    
else
    o.Loads.SectionID={[]};
    o.Loads.DeviceNumber={[]};
    o.Loads.LoadType={[]};
end
o.Loads.Phase1='X';
o.Loads.Phase2='Y';
o.Loads.Phase3='Z';

if isfield(o,'line') && isfield(o,'cable')
    % Process line and cable data to match LineCode structure created later on
    OHC=o.line;
    UGC=o.cable;
    % Combinining line codes for cables and lines; have to be careful with as
    % mixing cables and lines will likely result in non-unique ID (e.g.,
    % ID = 'default'), need to use the 'OH' and 'UG' line type designation to
    % sort this out later on (e.g., if 'OH' use first occurence, o/w use second
    % occurence of ID)
    o.LineCode.ConductorID=[OHC.LineID; UGC.CableID]; 
    o.LineCode.R1=[str2double(OHC.R1); str2double(UGC.R1)];
    o.LineCode.R0=[str2double(OHC.R0); str2double(UGC.R0)];
    o.LineCode.X1=[str2double(OHC.X1); str2double(UGC.X1)];
    o.LineCode.X0=[str2double(OHC.X0); str2double(UGC.X0)];
    o.LineCode.B1=[str2double(OHC.B1); str2double(UGC.B1)];
    o.LineCode.B0=[str2double(OHC.B0); str2double(UGC.B0)];
    o.LineCode.Z1=sqrt(o.LineCode.R1.^2+o.LineCode.X1.^2);
    o.LineCode.Z0=sqrt(o.LineCode.R0.^2+o.LineCode.X0.^2);
    o.LineCode.Rs=1/3*(2*o.LineCode.R1+o.LineCode.R0);
    o.LineCode.Xs=1/3*(2*o.LineCode.X1+o.LineCode.X0);
    o.LineCode.Zs=1/3*(2*o.LineCode.Z1+o.LineCode.Z0);
    o.LineCode.B1=[str2double(OHC.B1); str2double(UGC.B1)]; % uS/mile, currently not used 
    o.LineCode.B0=[str2double(OHC.B0); str2double(UGC.B0)];
    o.LineCode.R1=num2cell(o.LineCode.R1);
    o.LineCode.R0=num2cell(o.LineCode.R0);
    o.LineCode.X1=num2cell(o.LineCode.X1);
    o.LineCode.X0=num2cell(o.LineCode.X0);
    o.LineCode.B1=num2cell(o.LineCode.B1);
    o.LineCode.B0=num2cell(o.LineCode.B0);
    o.LineCode.Z1=num2cell(o.LineCode.Z1);
    o.LineCode.Z0=num2cell(o.LineCode.Z0);
end


% map line information to .Section
% it appears that line information comes in two flavors:
% (1)   line info is in .overheadlinesetting and .undergroundlinesetting,
%       in that case the code merges the info from the two objects to the
%       .Section and marks the line type as either 'OH', 'UG', or 'other'
%       not sure what 'other' signifies - perhaps a switch?
% (2)   line info is in .lineconfiguration and .lineconfiguration2, the
%       former has the line impedances and we can deal with it in a similar
%       fashion as the .overheadlinesetting and .undergroundlinesetting
%       objects, the latter specifies the lines through geometry, OpenDSS
%       can deal with that, but this will be tricky to convert to EMTP as
%       it requires calculating the impedances by using the auxiliary
%       program 'lineconstant' (or something similar)

% Merge .overheadlinesetting and .undergroundcablesetting, and put into .Section
if isfield(o,'overheadlinesetting') && isfield(o,'undergroundlinesetting')
    SID=o.Section.SectionID;
    OH=o.overheadlinesetting;
    UG=o.undergroundlinesetting;
    for i=1:length(SID)
        i_OH=GetIndex_StrInCell(OH.SectionID,SID{i});
        i_UG=GetIndex_StrInCell(UG.SectionID,SID{i});
        if length(i_OH)> 1 ||  length(i_UG)>1
            error('SectionID ' & SID{i} & ' has multiple underground or overhead lines assigned.')
        end
        if ~isempty(i_OH) && isempty(i_UG)
            o.Section.DeviceNumber{i,1}=OH.DeviceNumber{i_OH};
            o.Section.LineCableID{i,1}=OH.LineCableID{i_OH};
            o.Section.Length{i,1}=OH.Length{i_OH};
            o.Section.ConnectionStatus{i,1}=OH.ConnectionStatus{i_OH};
            o.Section.Phasing{i,1}=o.Section.Phasing{i};
            o.Section.LineType{i,1}='OH'; % added field for linetype info, OH=Overhead line, UG=Underground cable
        elseif isempty(i_OH) && ~isempty(i_UG)
            o.Section.DeviceNumber{i,1}=UG.DeviceNumber{i_UG};
            o.Section.LineCableID{i,1}=UG.LineCableID{i_UG};
            o.Section.Length{i,1}=UG.Length{i_UG};
            o.Section.ConnectionStatus{i,1}=UG.ConnectionStatus{i_UG};
            o.Section.Phasing{i,1}=o.Section.Phasing{i};
            o.Section.LineType{i,1}='UG'; % added field for linetype info, OH=Overhead line, UG=Underground cable
        elseif isempty(i_OH) && isempty(i_UG)
            o.Section.DeviceNumber{i,1}='';
            o.Section.LineCableID{i,1}='';
            o.Section.Length{i,1}='';
            o.Section.ConnectionStatus{i,1}='';
            o.Section.Phasing{i,1}=o.Section.Phasing{i};
            o.Section.LineType{i,1}='other';
        elseif ~isempty(i_OH) && ~isempty(i_UG)
            error('SectionID ' & SID{i} & ' has a duplicate assignment (underground and overhead line).')
        end
    end
end

if isfield(o,'lineconfiguration') && isfield(o,'lineconfiguration2')
    SID=o.Section.SectionID;
    OH=o.lineconfiguration;
    OH2=o.lineconfiguration2;
    for i=1:length(SID)
        i_OH=GetIndex_StrInCell(OH.SectionID,SID{i});
        i_OH2=GetIndex_StrInCell(OH2.SectionID,SID{i});
        if ~isempty(i_OH)
            if length(i_OH)> 1 
                error('SectionID ' & SID{i} & ' has multiple underground or overhead lines assigned.')
            end
            o.Section.DeviceNumber{i,1}=OH.DeviceNumber{i_OH};
            o.Section.CondID_A{i,1}=[];
            o.Section.CondID_B{i,1}=[];
            o.Section.CondID_C{i,1}=[];
            o.Section.CondID_N{i,1}=[];
            o.Section.LineCableID{i,1}=OH.LineCableID{i_OH};
            o.Section.Length{i,1}=OH.Length{i_OH};
            o.Section.ConnectionStatus{i,1}=OH.ConnectionStatus{i_OH};
            o.Section.SpacingID{i,1}=[];
            o.Section.Phasing{i,1}=o.Section.Phasing{i};
            if strcmp(OH.Overhead(i_OH),'1')
                o.Section.LineType{i,1}='OH'; % added field for linetype info, OH=Overhead line, UG=Underground cable
            else
                o.Section.LineType{i,1}='UG';
            end
        elseif ~isempty(i_OH2)
            o.Section.DeviceNumber{i,1}=OH2.DeviceNumber{i_OH2};
            o.Section.CondID_A{i,1}=OH2.CondID_A{i_OH2};
            o.Section.CondID_B{i,1}=OH2.CondID_B{i_OH2};
            o.Section.CondID_C{i,1}=OH2.CondID_C{i_OH2};
            o.Section.CondID_N{i,1}=OH2.CondID_N{i_OH2};
            o.Section.Length{i,1}=OH2.Length{i_OH2};
            o.Section.ConnectionStatus{i,1}=OH2.ConnectionStatus{i_OH2};
            o.Section.SpacingID{i,1}=OH2.SpacingID{i_OH2};
            o.Section.Phasing{i,1}=o.Section.Phasing{i};
            o.Section.LineType{i,1}='geometry';
        else
            o.Section.DeviceNumber{i,1}=[];
            o.Section.CondID_A{i,1}=[];
            o.Section.CondID_B{i,1}=[];
            o.Section.CondID_C{i,1}=[];
            o.Section.CondID_N{i,1}=[];
            o.Section.Length{i,1}=[];
            o.Section.ConnectionStatus{i,1}=[];
            o.Section.SpacingID{i,1}=[];
            o.Section.Phasing{i,1}=o.Section.Phasing{i};
            o.Section.LineType{i,1}='switch';
        end
    end
end


% Add information about shunt caps
if ~isempty(o.Capacitors)
    SID=o.Capacitors.SectionID;
    for i=1:length(SID)
        i_C=GetIndex_StrInCell(o.Section.SectionID,SID{i}); 
        o.Capacitors.bus1(i,1)=o.Section.FromNode(i_C);
        o.Capacitors.bus2(i,1)=o.Section.ToNode(i_C);
    end
end

% Add switch information
if ~isempty(o.Switches)
    sws=o.Switches;
    for i=1:length(sws.SectionID)
        if strcmpi(sws.ClosedPhase(i),'none')
            o.Switches.SwitchIsOpen(i,1)={1};
        else
            o.Switches.SwitchIsOpen(i,1)={0};
        end
    end
end

% Add fuse information
if ~isempty(o.Fuses)
    fs=o.Fuses;
    for i=1:length(fs.SectionID)
        if strcmpi(fs.ClosedPhase(i),'none')
            o.Fuses.FuseIsOpen(i,1)={1};
        else
            o.Fuses.FuseIsOpen(i,1)={0};
        end
        o.Fuses.UniqueDeviceID(i,1) = {['f_' fs.SectionID{i}]};
        o.Fuses.NearFromNode=1; % Assume all fuses are at the from node of the section, o/w enter '0'
    end
end


%% Convert structure
% Data is in cell array that is part of a structure, 
% Benandu code require different format, i.e., each data point is an
% instance in the structure
ob.Section=structconv(o.Section);
if ~isempty(o.Capacitors); ob.Capacitors=structconv(o.Capacitors); end
ob.Loads = structconv(o.Loads);
ob.customerclass = structconv(o.customerclass);
ob.feeder = structconv(o.feeder);
ob.sourceequivalent = structconv(o.sourceequivalent);
ob.Switches = structconv(o.Switches);
ob.breaker = structconv(o.breaker);
ob.recloser = structconv(o.recloser);
ob.Node = structconv(o.Node);
ob.source = structconv(o.source);
ob.Fuses = structconv(o.Fuses);
ob.Substation = structconv(o.substation);
ob.LineCode = structconv(o.LineCode);

if isfield(o,'Transformers')
    ob.Transformers = structconv(o.Transformers);
end

if isfield(o,'overheadlinesetting') && isfield(o,'undergroundlinesetting')
    ob.undergroundlinesetting = structconv(o.undergroundlinesetting);
    ob.overheadlinesetting = structconv(o.overheadlinesetting);
end


