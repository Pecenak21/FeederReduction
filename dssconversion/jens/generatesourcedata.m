function [struct] = generatesourcedata(o)
% generate structure that assembles/calculates all the source data we need


% data from setting structure (describes network)
SourceID_setting = o.SourceID;
NodeID=o.NodeID;
NetworkID=o.NetworkID;


% data from equipment structure (equipment library)

if ~strcmp(o.MVA,'')
    MVA=str2double(o.MVA);
    flgPU=1;
else
    flgPU=0; % No MVA specified=>assuming impedances are actual (not pu)
end
SourceID_equipment = o.ID;

KVLL=str2double(o.KVLL);
R1=str2double(o.R1);
X1=str2double(o.X1);
R0=str2double(o.R0);
X0=str2double(o.X0);
if isfield(o,'ByPhVoltDegPh1') && isfield(o,'ByPhVoltDegPh2') % information from SynerGEE file, not sure if CYME has this info
    ByPhVoltDegPh1=str2double(o.ByPhVoltDegPh1);
    ByPhVoltDegPh2=str2double(o.ByPhVoltDegPh2);
end
Conn=o.Conn;
% PhaseAngle=o.PhaseAngle;
% TxfoConnection=o.TxfoConnection;
% ImpedanceUnit=o.ImpedanceUnit;

% take the section id from the setting structure and fill in the blanks
% add setting info
struct.SourceID=SourceID_setting;
struct.NodeID=NodeID;
struct.NetworkID=NetworkID;

for i=1:length(SourceID_setting)
    % add equipment info
%     [flg index]=ismember(SourceID_setting(i),SourceID_equipment);
%     struct.MVA(i,1)=MVA(index); % for some strange reason, using this MVA base results in inconsistent results wrt to the SC currents provided by the customer (HOL Short Circuit Results_sorted.xlsx)
    if flgPU % impedances are given in pu
        struct.MVA(i,1)=MVA(i); % for the HydroOttawa feeder, 100 MVA gives consistent SC currents, i.e., we can match the SC currents provided by the customer
        struct.KVLL(i,1)=KVLL(i);
        struct.R1_pu(i,1)=R1(i); % for HydroOttawa feeder, customer provided different set of impedances, but we did not get the SC currents provided by the customer with the impedances provided by the customer
        struct.X1_pu(i,1)=X1(i); % the impedances from the [substation] CYME file seem to be the correct ones (if we use a 100 MVA base)
        struct.R0_pu(i,1)=R0(i);
        struct.X0_pu(i,1)=X0(i);
        % calculate missing parameters
        struct.Z_base(i,1)=(struct.KVLL(i,1))^2/struct.MVA(i,1); % some references use VLN, but this seems to be incorrect (or at least does not apply here)
        struct.R1(i,1)=struct.R1_pu(i,1)*struct.Z_base(i,1);
        struct.X1(i,1)=struct.X1_pu(i,1)*struct.Z_base(i,1);
        struct.R0(i,1)=struct.R0_pu(i,1)*struct.Z_base(i,1);
        struct.X0(i,1)=struct.X0_pu(i,1)*struct.Z_base(i,1);
        
         % fill in some parameters that are needed in the Excel sheet
        struct.R2_pu(i,1)=struct.R1_pu(i,1);
        struct.X2_pu(i,1)=struct.X1_pu(i,1);
    else
        struct.KVLL(i,1)=KVLL(i);
        struct.R1(i,1)=R1(i); 
        struct.X1(i,1)=X1(i); 
        struct.R0(i,1)=R0(i);
        struct.X0(i,1)=X0(i);
    end
    
    
    
    
    struct.Conn(i,1)=Conn(i);
%     struct.PhaseAngle(i,1)=PhaseAngle(i);
%     struct.TxfoConnection(i,1)=TxfoConnection(i);
%     struct.ImpedanceUnit(i,1)=ImpedanceUnit(i);
    
    if isfield(o,'ByPhVoltDegPh1') && isfield(o,'ByPhVoltDegPh2')
        if(mod(ByPhVoltDegPh1 - ByPhVoltDegPh2,360) > 0) % Phase 1 leads
            struct.Sequence = 'pos';
        else
            struct.Sequence = 'neg';
        end
    end


    
   
    % fill in some parameters that are needed in the Excel sheet
    struct.InServiceFlag(i,1)=1; % if it's listed, it is in service
    struct.R2(i,1)=struct.R1(i,1);
    struct.X2(i,1)=struct.X1(i,1);
    struct.RN(i,1)=0;
    struct.XN(i,1)=0;
    struct.PeakLoad(i,1)=0;
    struct.LightLoad(i,1)=0;
    

end

end



