function [ ob,node_bk ] = cyme2obj( o,flgTest,NetworkID)

% Convert cyme files with multiple sheets to object oriented struct
if nargin<2
    flgTest=0;
end


%% Convert structures content into cell arrays

% Load data
FieldNames_load=fieldnames(o.load);
for i=1:length(FieldNames_load)
    o.load.(FieldNames_load{i})=structconv(o.load.(FieldNames_load{i}));
    if isempty(o.load.(FieldNames_load{i}))
        o.load=rmfield(o.load,FieldNames_load{i});
    end
end

% Network data
o.network.feeder=o.network.section.feeder;
o.network.section=rmfield(o.network.section,'feeder');
o.network.section=o.network.section.section;
FieldNames_network=fieldnames(o.network);
for i=1:length(FieldNames_network)
    [r c]=size(o.network.(FieldNames_network{i}));
    if c==1
        o.network.(FieldNames_network{i})=structconv(o.network.(FieldNames_network{i}));
        if isempty(o.network.(FieldNames_network{i}))
            o.network=rmfield(o.network,FieldNames_network{i});
        end
    else
        for i_=1:c
            o.network.([FieldNames_network{i} num2str(i_)])=structconv(o.network.(FieldNames_network{i}){i_});
        end
        o.network=rmfield(o.network,FieldNames_network{i});
    end
end

% Equipment data
FieldNames_equipment=fieldnames(o.equipment);
for i=1:length(FieldNames_equipment)
    o.equipment.(FieldNames_equipment{i})=structconv(o.equipment.(FieldNames_equipment{i}));
    if isempty(o.equipment.(FieldNames_equipment{i}))
        o.equipment=rmfield(o.equipment,FieldNames_equipment{i});
    end
end


%% Reduce data

if flgTest
    % Use test system defined by Section IDs in the TestSystem.txt file
%     o.network.section.section.SectionID=o.TestSystem{:};
    [flg,index]=ismember(o.network.section.SectionID,o.TestSystem);
    index=find(index);
    o.network.section=structreduce(o.network.section,index);
end

% Select the system associated with the specified NetworkID
if ~strcmp(NetworkID,'')
    [flg,index]=ismember(o.network.section.NetworkID,NetworkID);
    index=find(index);
    o.network.section=structreduce(o.network.section,index);
    [flg,index]=ismember(o.network.headnodes.NetworkID,NetworkID);
    index=find(index);
    o.network.headnodes=structreduce(o.network.headnodes,index);
end

% only include objects and buses that are associated with a Section ID
% This is needed when dealing with smaller test sytems, but it is also good
% for getting rid of any "loose" objects and buses.

if ~isfield(o.network,'node') % create .node if multiple node types exist
    NodeID= [o.network.node1.NodeID; o.network.node2.NodeID];
    o.network.node.NodeID=NodeID;
    % in the SCE system, multiple nodes exist because some nodes have two x
    % coordinates and two y coordinates, here we make the nodes equal by
    % picking the first x coordinate and the first y coordinate.
    if isfield(o.network.node1,'CoordX')
        CoordXa=o.network.node1.CoordX;
    else
        CoordXa=o.network.node1.CoordX1;
    end
    if isfield(o.network.node2,'CoordX')
        CoordXb=o.network.node2.CoordX;
    else
        CoordXb=o.network.node2.CoordX1;
    end
    CoordX=[CoordXa; CoordXb];
    if isfield(o.network.node1,'CoordY')
        CoordYa=o.network.node1.CoordY;
    else
        CoordYa=o.network.node1.CoordY1;
    end
    if isfield(o.network.node2,'CoordY')
        CoordYb=o.network.node2.CoordY;
    else
        CoordYb=o.network.node2.CoordY1;
    end
    CoordY=[CoordYa; CoordYb];
    o.network.node.CoordX=CoordX;
    o.network.node.CoordY=CoordY;
end

% reduce nodes
[flg,index_FromNode]=ismember(o.network.section.FromNodeID,o.network.node.NodeID);
[flg,index_ToNode]=ismember(o.network.section.ToNodeID,o.network.node.NodeID);
index=[index_ToNode; index_FromNode];
% find the node the source connects to and make sure that that one is
% retained
[flg,index_SourceNode_FromNode]=ismember(o.network.source.NodeID,o.network.section.FromNodeID);
i_=find(index_SourceNode_FromNode);
if ~isempty(i_)
    index_SourceNode_section=index_SourceNode_FromNode(i_);
    [flg,index_SourceNode_node]=ismember(o.network.section.FromNodeID(index_SourceNode_section),o.network.node.NodeID);
    index=[index; index_SourceNode_node];
end
[flg,index_SourceNode_ToNode]=ismember(o.network.source.NodeID,o.network.section.ToNodeID);
i_=find(index_SourceNode_ToNode);
if ~isempty(i_)
    index_SourceNode_section=index_SourceNode_ToNode(i_);
    [flg,index_SourceNode_node]=ismember(o.network.section.ToNodeID(index_SourceNode_section),o.network.node.NodeID);
    index=[index; index_SourceNode_node];
end

node_bk = o.network.node;%VZ:stores the original nodes with IDs and Coordinates. needed later to add reduced sections.
o.network.node=structreduce(o.network.node,index);

% add the headnodes to the node structure, headnodes are where the
% substation sources are connected to
% o.network=structadd(o.network,'headnodes','NodeID',o.network,'node','NodeID','source','ZoneID');
% % the x and y coordinates of the headnodes reside in the .node objects,
% % pull this information from there and add
% [flg index]=ismember(o.network.node.CoordX,'');
% i_=find(index);
% 
% for i=1:length(i_)
%     k=i_(i);
%     [flg index]=ismember(o.network.node.NodeID,o.network.node.CoordX);
%     o.network.node.CoordX(i_(i))=;
% end
% o.network=rmfield(o.network,'headnodes');

% reduce content in .network
FieldNames=fieldnames(o.network);
for i=1:length(FieldNames)
    if isfield(o.network.(FieldNames{i}),'SectionID') % reduce objects with Section IDs
        [flg,index]=ismember(o.network.section.SectionID,o.network.(FieldNames{i}).SectionID);
%         index=find(index);
        o.network.(FieldNames{i})=structreduce(o.network.(FieldNames{i}),index);
        if isempty(o.network.(FieldNames{i}).SectionID)
            o.network=rmfield(o.network,FieldNames{i});
        end
    elseif isfield(o.network.(FieldNames{i}),'NodeID') % reduce objects with Node IDs
        [flg,index]=ismember(o.network.node.NodeID,o.network.(FieldNames{i}).NodeID);
% %         index=find(index);
        o.network.(FieldNames{i})=structreduce(o.network.(FieldNames{i}),index);
        if isempty(o.network.(FieldNames{i}).NodeID)
            o.network=rmfield(o.network,FieldNames{i});
        end
    elseif ~isfield(o.network.(FieldNames{i}),'ID') && ~isfield(o.network.(FieldNames{i}),'EqID')
        o.network=rmfield(o.network,FieldNames{i});
    end
end

% reduce content in .load, be careful not to reduce phases B and C away
FieldNames=fieldnames(o.load);
for i=1:length(FieldNames)
    if isfield(o.load.(FieldNames{i}),'SectionID') % reduce objects with Section IDs
        [flg,index]=ismember(o.load.(FieldNames{i}).SectionID,o.network.section.SectionID); % Careful! Non-unique section ID because one section can have many loads connected to it
        index=find(index);
        o.load.(FieldNames{i})=structreduce(o.load.(FieldNames{i}),index);
        if isempty(o.load.(FieldNames{i}).SectionID)
            o.load=rmfield(o.load,FieldNames{i});
        end
    elseif isfield(o.load.(FieldNames{i}),'NodeID') % reduce objects with Node IDs
        [flg,index]=ismember(o.network.node.NodeID,o.load.(FieldNames{i}).NodeID);
%         index=find(index);
        o.load.(FieldNames{i})=structreduce(o.load.(FieldNames{i}),index);
        if isempty(o.load.(FieldNames{i}).NodeID)
            o.load=rmfield(o.load,FieldNames{i});
        end
    elseif ~isfield(o.load.(FieldNames{i}),'ID') && ~isfield(o.load.(FieldNames{i}),'EqID')
        o.load=rmfield(o.load,FieldNames{i});
    end
end

%% NODES
o.Nodes=o.network.node;


%% TRANSFORMERS

if isfield(o.network,'transformer_setting') && isfield(o.equipment,'transformer')
    o.Transformers=structmerge(o.network.transformer_setting,'EqID',o.equipment.transformer,'ID',1,1);
    o.Transformers=structmerge(o.Transformers,'SectionID',o.network.section,'SectionID',1,1);

    % add processed info
%     o.Transformers.Phases=o.Transformers.Phase;
    o.Transformers.HighSideNearFromNode(1:length(o.Transformers.SectionID),1)=deal({'1'});
    connection = regexp(o.Transformers.Conn, '_', 'split');
    for i=1:length(connection) % js: This is rather clumsy, but I do not know another way to pull the information out of the cells that contain two cells
        o.Transformers.HighSideConnectionCode(i,1)=connection{i}(1);
        o.Transformers.LowSideConnectionCode(i,1)=connection{i}(2);
        o.Transformers.UniqueDeviceID(i,1) = {['t_' o.Transformers.SectionID{i}]};
    end
        
    % clean up
    o.network=rmfield(o.network,'transformer_setting');
    o.equipment=rmfield(o.equipment,'transformer');
% else
%     o.Transformers.SectionID='';
end

% add transformer information to Section Structure
o=structadd(o,'Transformers','SectionID',o,'Section','SectionID','transformer','SectionType');


%% Regulators

if isfield(o,'Regulators') % SynerGEE object
    o.Regulators=structmerge(o.Regulators,'SectionId',o.network.section,'SectionID',1,1);
% else
%     o.Regulators.SectionID='';
end

%% LOADS
% Add information from .customer_class and .customer_loads to .LoadsAll
if isfield(o.load,'customer_class') && isfield(o.load,'customer_loads')
    o.LoadsAll=structmerge(o.load.customer_loads,'CustomerType',o.load.customer_class,'ID',1,1);
    o.load=rmfield(o.load,'customer_class');
    o.load=rmfield(o.load,'customer_loads');
end

% Add information from .large_customer_loads to .LoadsAll
if isfield(o.load,'large_customer_loads')
    fn_LA=fieldnames(o.LoadsAll);
    fn_lc=fieldnames(o.load.large_customer_loads);
    n=length(o.LoadsAll.SectionID);
    n2=length(o.load.large_customer_loads.SectionID);
    for i=1:length(fn_lc)
        if ismember(fn_lc(i),fn_LA)
            o.LoadsAll.(fn_lc{i})=[o.LoadsAll.(fn_lc{i});o.load.large_customer_loads.(fn_lc{i})];
        else
            o.LoadsAll.(fn_lc{i})=[cell(n,1);o.load.large_customer_loads.(fn_lc{i})];
        end
    end
    for i=1:length(fn_LA)
        if ~ismember(fn_LA(i),fn_lc)
            o.LoadsAll.(fn_LA{i})=[o.LoadsAll.(fn_LA{i});cell(n2,1)];
        end
    end
end

% % Add information from .loads to .Loads
% if isfield(o.load,'loads') && isfield(o.Loads,'SectionID')
%     o.Loads=structmerge(o.Loads,'SectionID',o.load.loads,'SectionID',1,1);
%     o.load=rmfield(o.load,'loads');
% end

% add processed info
if isfield(o,'LoadsAll')
    SID=unique(o.LoadsAll.SectionID);
    for i=1:length(SID)
        o.Loads.SectionID{i,1}=SID{i};
        if isfield(o.LoadsAll,'LocationToModelSpotLoads')
        	o.Loads.LocationToModelSpotLoads{i,1}=o.LoadsAll.LocationToModelSpotLoads{i,1}; % SynerGEE specifies this and the info is used here, CYME does not seem to specify this
        else
        	o.Loads.LocationToModelSpotLoads{i,1}='C'; % default is that load is at 'to' node, change to 'F' if load is at from node
        end
%         o.Loads.LocationToModelSpotLoads{i,1}='F'; % 'F' = load is modeled at the 'from' node, otherwise load is modeled at the 'to' node
        if strcmpi(o.LoadsAll.LoadType(i),'spot')
        	o.Loads.IsSpotLoad{i,1}='1';
        else
        	o.Loads.IsSpotLoad{i,1}='0';
        end
        o.Loads.Phase1{i,1}='';
        o.Loads.Phase2{i,1}='';
        o.Loads.Phase3{i,1}='';
        o.Loads.PhaseCheck{i,1}='';
        [flg index]=ismember(o.LoadsAll.SectionID,SID(i));
        index=find(index); % this gives the indeces of all loads associated with SID(i)
        o.Loads.ConstantPower{i,1}=o.LoadsAll.ConstantPower{index(1)};
        o.Loads.ConstantCurrent{i,1}=o.LoadsAll.ConstantCurrent{index(1)};
        o.Loads.ConstantImpedance{i,1}=o.LoadsAll.ConstantImpedance{index(1)};
        
        % SynerGEE stuff
        if isfield(o.LoadsAll,'GenTotalKw')
            o.Loads.GenTotalKw{i,1}=o.LoadsAll.GenTotalKw{index(1)};
            o.Loads.GenTotalKvar{i,1}=o.LoadsAll.GenTotalKvar{index(1)};
            o.Loads.GenStatus{i,1}=o.LoadsAll.GenStatus{index(1)};
            o.Loads.GenPct{i,1}=o.LoadsAll.GenPct{index(1)};
            o.Loads.GenCustClass{i,1}=o.LoadsAll.GenCustClass{index(1)};
            o.Loads.UniqueDeviceId{i,1}=o.LoadsAll.UniqueDeviceId{index(1)};
        end
        
        % Phase ABC - sometimes the loads come as ABC. These need to be included in the list as well
        [flg2 index2]=ismember(o.LoadsAll.LoadPhase(index),'ABC');
        if (flg2~=0)
            index2=index.*index2;
            index2=index2(index2~=0);
            PhaseKVA=sum(str2double(o.LoadsAll.Value1(index2)));
            PF=str2double(o.LoadsAll.Value2(index2))/100;
            PhaseKw=PhaseKVA*PF;
            PhaseKvar=sqrt(PhaseKw^2/PF^2-PhaseKw^2);
            o.Loads.Phase1Kw{i,1}=num2str(PhaseKw/3);o.Loads.Phase2Kw{i,1}=num2str(PhaseKw/3);o.Loads.Phase3Kw{i,1}=num2str(PhaseKw/3);
            o.Loads.Phase1Kvar{i,1}=num2str(PhaseKvar/3);o.Loads.Phase2Kvar{i,1}=num2str(PhaseKvar/3);o.Loads.Phase3Kvar{i,1}=num2str(PhaseKvar/3);
            o.Loads.Phase1Kwh{i,1}=num2str(str2double(o.LoadsAll.KWH(index2))/3);
            o.Loads.Phase1ConnectedKVA{i,1}=num2str(str2double(o.LoadsAll.ConnectedKVA(index2))/3);
            o.Loads.Phase1Customers{i,1}=num2str(str2double(o.LoadsAll.NumberOfCustomer(index2))/3);
            o.Loads.Phase1{i,1}='X';
            o.Loads.Phase2Kwh{i,1}=num2str(str2double(o.LoadsAll.KWH(index2))/3);
            o.Loads.Phase2ConnectedKVA{i,1}=num2str(str2double(o.LoadsAll.ConnectedKVA(index2))/3);
            o.Loads.Phase2Customers{i,1}=num2str(str2double(o.LoadsAll.NumberOfCustomer(index2))/3);
            o.Loads.Phase2{i,1}='Y';
            o.Loads.Phase3Kwh{i,1}=num2str(str2double(o.LoadsAll.KWH(index2))/3);
            o.Loads.Phase3ConnectedKVA{i,1}=num2str(str2double(o.LoadsAll.ConnectedKVA(index2))/3);
            o.Loads.Phase3Customers{i,1}=num2str(str2double(o.LoadsAll.NumberOfCustomer(index2))/3);
            o.Loads.Phase3{i,1}='Z';
            o.Loads.PhaseCheck{i,1}='XYZ';
            continue
        end
        
        % Phase A
        [flg2 index2]=ismember(o.LoadsAll.LoadPhase(index),'A');
        index2=index.*index2;
        index2=index2(index2~=0); % this gives the indeces of all phase A loads for SID(i)
        if ~isempty(index2)
%             Phase1Kw=sum(str2double(o.LoadsAll.Value1(index2)));
%             o.Loads.Phase1Kw{i,1}=num2str(Phase1Kw);
            Phase1KVA=str2double(o.LoadsAll.Value1(index2));%VZ: Sunil confirmed that value 1 is kVA not kW
            Phase1KVA=Phase1KVA(Phase1KVA~=0);
            Value2=str2double(o.LoadsAll.Value2(index2));
            Value2=Value2(Value2~=0);
            if ~isempty(Value2)
                Phase1Value2=str2double(o.LoadsAll.Value2(index2));
                Phase1Value2=Phase1Value2(Phase1Value2~=0);
%                 n1Value2=length(o.LoadsAll.Value2(index2));
%                 n1Value2=nnz(str2double(o.LoadsAll.Value2(index2)));%VZ: count nonzero elements rather than all elements. sometimes same node and phase will have two loads one with kw=0&pf=0 and the other with some legal values. averaging both of them will falsify the outcome. 
                pf1=Phase1Value2/100; % pf=P/S; S^2=P^2+Q^2; Q=sqrt(P^2/pf^2-P^2)                
                if length(Value2)>1
                    Phase1Kw = 0;Phase1KVAr = 0;
                    for zz=1:length(Value2)
                       Kw = Phase1KVA(zz)*pf1(zz); 
                       KVAr = sqrt(Kw^2/pf1(zz)^2-Kw^2);
                       Phase1Kw = Phase1Kw + Kw;
                       Phase1KVAr = Phase1KVAr + KVAr;
                    end
                else
                    Phase1Kw = Phase1KVA*pf1;
                    Phase1KVAr = sqrt(Phase1Kw^2/pf1^2-Phase1Kw^2);
                end
                o.Loads.Phase1Kw{i,1}=num2str(Phase1Kw);
                o.Loads.Phase1Kvar{i,1}=num2str(Phase1KVAr);
            else
                o.Loads.Phase1Kvar{i,1}='0';
                o.Loads.Phase1Kw{i,1}='0';
            end
            o.Loads.Phase1Kwh{i,1}=num2str(sum(str2double(o.LoadsAll.KWH(index2))));
            o.Loads.Phase1ConnectedKVA{i,1}=num2str(sum(str2double(o.LoadsAll.ConnectedKVA(index2))));
            o.Loads.Phase1Customers{i,1}=num2str(sum(str2double(o.LoadsAll.NumberOfCustomer(index2))));
            if ~strcmp(o.Loads.Phase1Kw{i,1},'0') || ~strcmp(o.Loads.Phase1Kvar{i,1},'0')
                o.Loads.Phase1{i,1}='X';
                o.Loads.PhaseCheck{i,1}=[o.Loads.PhaseCheck{i,1} 'X'];
            end
        else
            o.Loads.Phase1Kw{i,1}='';
            o.Loads.Phase1Kvar{i,1}='';
            o.Loads.Phase1Kwh{i,1}='';
            o.Loads.Phase1ConnectedKVA{i,1}='';
            o.Loads.Phase1Customers{i,1}='';
            
        end
        
        % Phase B
        [flg2 index2]=ismember(o.LoadsAll.LoadPhase(index),'B');
        index2=index.*index2;
        index2=index2(index2~=0); % this gives the indeces of all phase B loads for SID(i)
        if ~isempty(index2)
%             Phase2Kw=sum(str2double(o.LoadsAll.Value1(index2)));
%             o.Loads.Phase2Kw{i,1}=num2str(Phase2Kw);
            Phase2KVA=str2double(o.LoadsAll.Value1(index2));%VZ: Sunil confirmed that value 1 is kVA not kW
            Phase2KVA=Phase2KVA(Phase2KVA~=0);
            Value2=str2double(o.LoadsAll.Value2(index2));
            Value2=Value2(Value2~=0);
            if ~isempty(Value2)
                Phase2Value2=str2double(o.LoadsAll.Value2(index2));
                Phase2Value2=Phase2Value2(Phase2Value2~=0);
%                 n2Value2=length(o.LoadsAll.Value2(index2));
%                 n2Value2=nnz(str2double(o.LoadsAll.Value2(index2)));%VZ: count nonzero elements rather than all elements. sometimes same node and phase will have two loads one with kw=0&pf=0 and the other with some legal values. averaging both of them will falsify the outcome.
                pf2=Phase2Value2/100; % pf=P/S; S^2=P^2+Q^2; Q=sqrt(P^2/pf^2-P^2)
                if length(Value2)>1
                    Phase2Kw = 0;Phase2KVAr = 0;
                    for zz=1:length(Value2)
                       Kw = Phase2KVA(zz)*pf2(zz); 
                       KVAr = sqrt(Kw^2/pf2(zz)^2-Kw^2);
                       Phase2Kw = Phase2Kw + Kw;
                       Phase2KVAr = Phase2KVAr + KVAr;
                    end
                else
                    Phase2Kw = Phase2KVA*pf2;
                    Phase2KVAr = sqrt(Phase2Kw^2/pf2^2-Phase2Kw^2);
                end
                o.Loads.Phase2Kw{i,1}=num2str(Phase2Kw);
                o.Loads.Phase2Kvar{i,1}=num2str(Phase2KVAr);
            else
                o.Loads.Phase2Kvar{i,1}='0';
                o.Loads.Phase2Kw{i,1}='0';
            end
            o.Loads.Phase2Kwh{i,1}=num2str(sum(str2double(o.LoadsAll.KWH(index2))));
            o.Loads.Phase2ConnectedKVA{i,1}=num2str(sum(str2double(o.LoadsAll.ConnectedKVA(index2))));
            o.Loads.Phase2Customers{i,1}=num2str(sum(str2double(o.LoadsAll.NumberOfCustomer(index2))));
            if ~strcmp(o.Loads.Phase2Kw{i,1},'0') || ~strcmp(o.Loads.Phase2Kvar{i,1},'0')
                o.Loads.Phase2{i,1}='Y';
                o.Loads.PhaseCheck{i,1}=[o.Loads.PhaseCheck{i,1} 'Y'];
            end
         else
            o.Loads.Phase2Kw{i,1}='';
            o.Loads.Phase2Kvar{i,1}='';
            o.Loads.Phase2Kwh{i,1}='';
            o.Loads.Phase2ConnectedKVA{i,1}='';
            o.Loads.Phase2Customers{i,1}='';
            
        end
        
        % Phase C
        [flg2 index2]=ismember(o.LoadsAll.LoadPhase(index),'C');
        index2=index.*index2;
        index2=index2(index2~=0); % this gives the indeces of all phase C loads for SID(i)
        if ~isempty(index2)
%             Phase3Kw=sum(str2double(o.LoadsAll.Value1(index2)));
%             o.Loads.Phase3Kw{i,1}=num2str(Phase3Kw);
            Phase3KVA=str2double(o.LoadsAll.Value1(index2));%VZ: Sunil confirmed that value 1 is kVA not kW
            Phase3KVA=Phase3KVA(Phase3KVA~=0);
            Value2=str2double(o.LoadsAll.Value2(index2));
            Value2=Value2(Value2~=0);
            if ~isempty(Value2)
                Phase3Value2=str2double(o.LoadsAll.Value2(index2));
                Phase3Value2=Phase3Value2(Phase3Value2~=0);
%                 n3Value2=length(o.LoadsAll.Value2(index2));
%                 n3Value2=nnz(str2double(o.LoadsAll.Value2(index2)));%VZ: count nonzero elements rather than all elements. sometimes same node and phase will have two loads one with kw=0&pf=0 and the other with some legal values. averaging both of them will falsify the outcome.
                pf3=Phase3Value2/100; % pf=P/S; S^2=P^2+Q^2; Q=sqrt(P^2/pf^2-P^2)
                if length(Value2)>1
                    Phase3Kw = 0;Phase3KVAr = 0;
                    for zz=1:length(Value2)
                       Kw = Phase3KVA(zz)*pf3(zz); 
                       KVAr = sqrt(Kw^2/pf3(zz)^2-Kw^2);
                       Phase3Kw = Phase3Kw + Kw;
                       Phase3KVAr = Phase3KVAr + KVAr;
                    end
                else
                    Phase3Kw = Phase3KVA*pf3;
                    Phase3KVAr = sqrt(Phase3Kw^2/pf3^2-Phase3Kw^2);
                end
                o.Loads.Phase3Kw{i,1}=num2str(Phase3Kw);
                o.Loads.Phase3Kvar{i,1}=num2str(Phase3KVAr);
            else
                o.Loads.Phase3Kvar{i,1}='0';
                o.Loads.Phase3Kw{i,1}='0';
            end
            o.Loads.Phase3Kwh{i,1}=num2str(sum(str2double(o.LoadsAll.KWH(index2))));
            o.Loads.Phase3ConnectedKVA{i,1}=num2str(sum(str2double(o.LoadsAll.ConnectedKVA(index2))));
            o.Loads.Phase3Customers{i,1}=num2str(sum(str2double(o.LoadsAll.NumberOfCustomer(index2))));
            if ~strcmp(o.Loads.Phase3Kw{i,1},'0') || ~strcmp(o.Loads.Phase3Kvar{i,1},'0')
                o.Loads.Phase3{i,1}='Z';
                o.Loads.PhaseCheck{i,1}=[o.Loads.PhaseCheck{i,1} 'Z'];
            end
        else
            o.Loads.Phase3Kw{i,1}='';
            o.Loads.Phase3Kvar{i,1}='';
            o.Loads.Phase3Kwh{i,1}='';
            o.Loads.Phase3ConnectedKVA{i,1}='';
            o.Loads.Phase3Customers{i,1}='';
            
        end
    end
    [flg,index]=ismember(o.Loads.PhaseCheck,'');
    index=find(~index);
    o.Loads=structreduce(o.Loads,index);
    o.Loads=structmerge(o.Loads,'SectionID',o.network.section,'SectionID',1,1);

end

%% LINES
% it appears that line information comes in two flavors:
% (1)   line info is in .overheadlinesetting and .undergroundlinesetting,
%       in that case the code merges the info from the two objects to the
%       .Section and marks the line type as either 'OH', 'UG', or 'other'
%       not sure what 'other' signifies - perhaps a switch?
%       This applies to the HydroOttawa system.
% (2)   line info is in .lineconfiguration and .lineconfiguration2, the
%       former has the line impedances and we can deal with it in a similar
%       fashion as the .overheadlinesetting and .undergroundlinesetting
%       objects, the latter specifies the lines through geometry, OpenDSS
%       can deal with that, but this will be tricky to convert to EMTP as
%       it requires calculating the impedances by using the auxiliary
%       program 'lineconstant' (or something similar)
%       This applies to the SCE system.

% get field names that are common to all objects that have line info
fn=GetFieldNames('intersect',o.equipment,'line',o.equipment,'cable',o.equipment,'concentric_neutral_cable');

% add line info to LineCode structure
for i=1:length(fn)
    o.LineCode.(fn{i})=[];
    if isfield(o.equipment,'line') 
        o.LineCode.(fn{i})=[o.LineCode.(fn{i}); o.equipment.line.(fn{i})];
    end
    if isfield(o.equipment,'cable') 
        o.LineCode.(fn{i})=[o.LineCode.(fn{i}); o.equipment.cable.(fn{i})];
    end
    if isfield(o.equipment,'concentric_neutral_cable') 
        o.LineCode.(fn{i})=[o.LineCode.(fn{i}); o.equipment.concentric_neutral_cable.(fn{i})];
    end
end

% add line geometry to structure
if isfield(o.equipment,'spacing_table_for_line')
    o.LineGeometry=o.equipment.spacing_table_for_line;
end

% add wire properties to structure
if isfield(o.equipment,'conductor')
    o.WireData=o.equipment.conductor;
end


% in HydroOttawa system the network information of the lines resides in the
% 'overheadline_setting' and 'undergroundline_setting' paragraphs in the
% network file
if isfield(o.network,'overheadline_setting') || isfield(o.network,'undergroundline_setting')
    % get field names that are common to all objects that have line info
    fn=GetFieldNames('intersect',o.network,'overheadline_setting',o.network,'undergroundline_setting');
    for i=1:length(fn)
        o.Line.(fn{i})=[];
        if isfield(o.network,'overheadline_setting') 
            o.Line_oh.(fn{i})=[o.Line.(fn{i}); o.network.overheadline_setting.(fn{i})];
        end
        if isfield(o.network,'undergroundline_setting') 
            o.Line_ug.(fn{i})=[o.Line.(fn{i}); o.network.undergroundline_setting.(fn{i})];
        end
    end
    if isfield(o,'Line_ug')
        % get rid of the Section IDs in .Line_ug that are already included
        % as transformers
        [flg index]=ismember(o.Line_ug.SectionID,o.Section.SectionID);
        index=find(~index);
        o.Line_ug=structreduce(o.Line_ug,index);
        % add additional content, this should not result in duplicate as I
        % trimmed the SectionIDs previously
        o=structadd(o,'Line_ug','SectionID',o,'Section','SectionID','UG','SectionType');
    end
    if isfield(o,'Line_oh')
        % get rid of the Section IDs in .Line_oh that are already included
        % as transformers
        [flg index]=ismember(o.Line_oh.SectionID,o.Section.SectionID);
        index=find(~index);
        o.Line_oh=structreduce(o.Line_oh,index);
        % add additional content, this should not result in duplicate as I
        % trimmed the SectionIDs previously
        o=structadd(o,'Line_oh','SectionID',o,'Section','SectionID','OH','SectionType');
        
    end

end

% in SCE system the network information of the lines resides in the
% 'line_configuration' paragraph in the network file
if isfield(o.network,'line_configuration1') || isfield(o.network,'line_configuration2') || isfield(o.network,'line_configuration')
    if isfield(o.network,'line_configuration')
        % add lines that are specified through SectionID in
        % line_configuration
        o=structadd(o.network,'line_configuration','SectionID',o,'Section','SectionID','OHUG','SectionType');
        clear index;
        index=ismember(o.Section.SectionType,{'OHUG'});
        o.Section.SectionType(index)=o.Section.Overhead(index);
        o.Section.SectionType=ReplaceStringsInCellArray(o.Section.SectionType,'1','OH');
        o.Section.SectionType=ReplaceStringsInCellArray(o.Section.SectionType,'0','UG');
    end
    
    if isfield(o.network,'line_configuration1')
        % add lines that are specified through SectionID in line_configuration1
        o=structadd(o.network,'line_configuration1','SectionID',o,'Section','SectionID','OHUG','SectionType');
        o.Section.SectionType=o.Section.Overhead;
        o.Section.SectionType=ReplaceStringsInCellArray(o.Section.SectionType,'1','OH');
        o.Section.SectionType=ReplaceStringsInCellArray(o.Section.SectionType,'0','UG');
        if isfield(o.network,'line_configuration2')
            % add lines that are specified through SectionID in line_configuration2
            o=structadd(o.network,'line_configuration2','SectionID',o,'Section','SectionID','geometry','SectionType');
        end
    end
end
    
% add loads that are specified through SectionID in Loads
% o=structadd(o,'Loads','SectionID',o,'Section','SectionID','load','SectionType');

% add default values for neutral grounding
o.Section.NeutIsGrounded(1:length(o.Section.SectionType),1)=deal({'1'});

%% SWITCHES
if isfield(o.network,'switch_setting')
%     o.Switches= o.network.switch_setting; % structadd(o.network,'switch_setting','SectionID',o,'Section','SectionID','switch','SectionType');
    fn=fieldnames(o.network.switch_setting);
    for i=1:length(fn)
        o.Switches.(fn{i})=o.network.switch_setting.(fn{i});
    end
    o.Switches=structmerge(o.Switches,'SectionID',o.network.section,'SectionID',1,1);
    % add switch information to Section Structure
    o=structadd(o,'Switches','SectionID',o,'Section','SectionID','switch','SectionType');
end

%% RECLOSERS
if isfield(o.network,'recloser_setting')
    o.Reclosers= o.network.recloser_setting; 
end


%% SHUNT CAPACITORS

% pull information from various structs
if isfield(o.network,'shunt_capacitor_setting')
    fn=fieldnames(o.network.shunt_capacitor_setting);
    for i=1:length(fn)
        o.Capacitors.(fn{i})=o.network.shunt_capacitor_setting.(fn{i});
    end
    o.Capacitors=structmerge(o.Capacitors,'SectionID',o.network.section,'SectionID',1,1);
    if isfield(o,'equipment')
        if isfield(o.equipment,'shunt_capacitor')
            o.Capacitors=structmerge(o.Capacitors,'ShuntCapacitorID',o.equipment.shunt_capacitor,'ID',1,0); % information in network struct is given preference
        end
    end
    % add cap information to Section Structure
    o=structadd(o,'Capacitors','SectionID',o,'Section','SectionID','capacitor','SectionType');
end



%% Source

% pull information from various structs
if isfield(o.network,'source')
    fn=fieldnames(o.network.source);
    for i=1:length(fn)
        o.Sources.(fn{i})=o.network.source.(fn{i});
    end
    o.Sources=structmerge(o.Sources,'SourceID',o.equipment.substation,'ID',1,0); % information in network struct is given preference
end




%% FUSES
% Section IDs for fuses are also section ID for lines, hence it should be
% OK to ignore fuses if not worried about protection stuff

% pull information from various structs
if isfield(o.network,'fuse_setting')
    fn=fieldnames(o.network.fuse_setting);
    for i=1:length(fn)
        o.Fuses.(fn{i})=o.network.fuse_setting.(fn{i});
    end
    o.Fuses=structmerge(o.Fuses,'SectionID',o.network.section,'SectionID',1,1);
    o.Fuses=structmerge(o.Fuses,'EqID',o.equipment.fuse,'ID',1,0); % information in network struct is given preference
end
% 
% %% BREAKERS
% % Section IDs for breakers are also section ID for lines, hence it should be
% % OK to ignore fuses if not worried about protection stuff
% 
% % pull information from various structs
% if isfield(o.network,'breaker_setting')
%     fn=fieldnames(o.network.breaker_setting);
%     for i=1:length(fn)
%         o.Breakers.(fn{i})=o.network.breaker_setting.(fn{i});
%     end
%     o.Breakers=structmerge(o.Breakers,'SectionID',o.network.section,'SectionID',1,1);
%     o.Breakers=structmerge(o.Breakers,'EqID',o.equipment.fuse,'ID',1,0); % information in network struct is given preference
% end

% %% SECTION
% % pull information from various structs
% if isfield(o.network,'section')
%     fn=fieldnames(o.network.section);
%     for i=1:length(fn)
%         o.Section.(fn{i})=o.network.section.(fn{i});
%     end
%     o.Breakers=structmerge(o.Breakers,'SectionID',o.network.section,'SectionID',1,1);
%     o.Breakers=structmerge(o.Breakers,'EqID',o.equipment.fuse,'ID',1,0); % information in network struct is given preference
% end

%% Clean up

% This gets rid of the empty fields and pulls NodeIDs from the original section structure
% (1) remove the NodeIDs because they have been only sporadicly added previosly
if isfield(o.Section,'FromNodeID');
    o.Section=rmfield(o.Section,'FromNodeID');
end
if isfield(o.Section,'ToNodeID');
    o.Section=rmfield(o.Section,'ToNodeID');
end
% (2) Pulls the NodeIDs back in from the o.network.section structure
o.Section=structmerge(o.Section,'SectionID',o.network.section,'SectionID',1,1);

% add mystery sections, i.e., sections that are neither associated with a
% line, switch, cap, fuse, nor transformer
% these section types are likely loads, which are dealt with separately
if isfield(o,'Transformers')
    index=~ismember(o.network.section.SectionID,[o.Section.SectionID;o.Transformers.SectionID]); % system may or may not have transformers
else
    index=~ismember(o.network.section.SectionID,o.Section.SectionID);
end
if any(index)
    fn=fieldnames(o.Section);
    MysterySections=o.network.section.SectionID(index);
    MysterySections_FromNodeID=o.network.section.FromNodeID(index);
    MysterySections_ToNodeID=o.network.section.ToNodeID(index);
    n_MysterySections=length(MysterySections);
    clear MysterySections_ClosedPhase;
    MysterySections_ClosedPhase(1:n_MysterySections,1)=deal({'ABC'});
%     o.Section.ClosedPhase=[o.Section.ClosedPhase;MysterySections_ClosedPhase];
    for i=1:length(fn)
        if strcmp(fn{i},'SectionID')
            o.Section.SectionID=[o.Section.SectionID;MysterySections];
        elseif strcmp(fn{i},'FromNodeID')
            o.Section.FromNodeID=[o.Section.FromNodeID;MysterySections_FromNodeID];
        elseif strcmp(fn{i},'ToNodeID')
            o.Section.ToNodeID=[o.Section.ToNodeID;MysterySections_ToNodeID];
        elseif strcmp(fn{i},'ClosedPhase')
            o.Section.ClosedPhase=[o.Section.ClosedPhase;MysterySections_ClosedPhase];
        elseif ~strcmp(fn{i},'SectionType')
            clear dummy;
            dummy(1:n_MysterySections,1)=deal({''});
            o.Section.(fn{i})=[o.Section.(fn{i});dummy];
        end
    end
    
    clear SectionType;
    SectionType(1:n_MysterySections,1)=deal({'other'});
    o.Section.SectionType=[o.Section.SectionType;SectionType];
    % implement mystery sections as three phase closed switches 

    
end

% get rid of switches
% o.Section=SectionTrim(o.Section,'SectionType','switch','FromNodeID','ToNodeID','FromNodeID');





%% Convert structure
% Data is in cell array that is part of a structure, 
% Benandu code require different format, i.e., each data point is an
% instance in the structure
ob.Section=structconv(o.Section);
if isfield(o,'Loads')
    ob.Loads=structconv(o.Loads);
end
if isfield(o,'LineCode')
    ob.LineCode=structconv(o.LineCode);
end
if isfield(o,'LineGeometry')
    ob.LineGeometry=structconv(o.LineGeometry);
end
if isfield(o,'WireData')
    ob.WireData=structconv(o.WireData);
end
if isfield(o,'Transformers')
    ob.Transformers = structconv(o.Transformers);
end
if isfield(o,'Regulators')
    ob.Regulators = structconv(o.Regulators);
end
if isfield(o,'Sources')
    ob.Sources = structconv(o.Sources);
end
if isfield(o,'Nodes')
    ob.Nodes = structconv(o.Nodes);
end
if isfield(o,'Fuses')
    ob.Fuses = structconv(o.Fuses);
end
if isfield(o,'Switches')
    ob.Switches = structconv(o.Switches);
end

if isfield(o,'Reclosers')
    ob.Reclosers = structconv(o.Reclosers);
end

if isfield(o,'Capacitors')
    ob.Capacitors = structconv(o.Capacitors);
end

% if ~isempty(o.shunt_capacitor_settings); ob.shunt_capacitor_settings=structconv(o.shunt_capacitor_settings); end
% ob.loads = structconv(o.loads);
% ob.customer_class = structconv(o.customer_class);
% ob.feeder = structconv(o.feeder);
% ob.sourceequivalent = structconv(o.sourceequivalent);
% ob.Switches = structconv(o.Switches);
% ob.breaker = structconv(o.breaker);
% ob.recloser_setting_setting = structconv(o.recloser_setting_setting);
% ob.node = structconv(o.node);
% ob.source = structconv(o.source);
% ob.fuse_settings = structconv(o.fuse_settings);
% ob.Substation = structconv(o.substation);
% ob.LineCode = structconv(o.LineCode);
% 
% if isfield(o,'Transformers')
%     ob.Transformers = structconv(o.Transformers);
% end
% 
% if isfield(o,'overheadlinesetting') && isfield(o,'undergroundlinesetting')
%     ob.undergroundlinesetting = structconv(o.undergroundlinesetting);
%     ob.overheadlinesetting = structconv(o.overheadlinesetting);
% end


