function[ExecutionFlag] = GetCYMData_network(FilePath)


% GetLinesFromFile(FilePath)
%
% PURPOSE : get content of each line in file
%
%
% INPUT :   FilePath: directory and file name
%           strFormat: format of data in file, 
%               e.g. '%s %f %f %f' 1st column is string , 2nd through 4th
%               columns are floating point numbers
%           strDelimiter: delimiter that separates columns in data file
%           flgHeader: 0 -> no header, 1 -> header
%           strDelimiter_Header: delimiter that separates header columns 
%
% OUTPUT :  y_First: First column in cell array
%           y: Complete cell array
%           y_Header: Header array (if flgHeader is set)
%           struc: a structure with header names and respective column content under each header name (if flgHeader is set)


yNODE=[];
yHEADNODES=[];
ySOURCE=[];
ySOURCEEQUIVALENT=[];
yLOADEQUIVALENT=[];
yOVERHEADLINESETTING=[];
yUNDERGROUNDLINESETTING=[];
ySECTION=[];
yTRANSFORMERSETTING=[];
ySWITCHSETTING=[];
yBREAKERSETTING=[];
yFUSESETTING=[];
yRECLOSERSETTING=[];
% ySECTIONALIZERSETTING=[];
% ySHUNTCAPACITORSETTING=[]; % HydroOttawa system does not appear to have caps
yINTERMEDIATENODES=[];
yDEVICETAG=[];
yAUTOTAPCHANGINGEXTENSION=[];


FileName=RemoveFileExtension(FilePath);

n_skip=9;

if ~exist(FilePath)
    str=[FilePath ' does not exist.'];
    msgbox(str)
    y_First='';
    y='';
    y_Header='';
    return
end

fid = fopen(FilePath, 'r');
i=1;
for i=1:n_skip % skip a n_skip lines
    tline = fgetl(fid);
end

% NODE data
fid2 = fopen([FileName '_NODE.txt'], 'w');
tline = fgetl(fid);
while isempty(strfind(tline,'[HEADNODES]'))
      yNODE=tline;
      fprintf(fid2,'%s \r\n',yNODE);
      tline = fgetl(fid);
end
fclose(fid2);

% HEADNODES data
fid2 = fopen([FileName '_HEADNODES.txt'], 'w');
if ~isempty(strfind(tline,'[HEADNODES]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_HEADNODES=NodeID,NetworkID'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SOURCE]')) 
         yHEADNODES=tline;
         fprintf(fid2,'%s\n',yHEADNODES);
         tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end
 
% Source data
fid2 = fopen([FileName '_SOURCE.txt'], 'w');
if ~isempty(strfind(tline,'[SOURCE]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SOURCE=SourceID,DeviceNumber,NodeID,ConnectorIndex,NetworkID,DesiredVoltage,StructureID,HarmonicEnveloppe,DemandType,Value1A,Value2A,Value1B,Value2B,Value1C,Value2C,Value1ABC,Value2ABC'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SOURCE EQUIVALENT]')) 
             ySOURCE=tline;
             fprintf(fid2,'%s\n',ySOURCE);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% Source data
fid2 = fopen([FileName '_SOURCE EQUIVALENT.txt'], 'w');
if ~isempty(strfind(tline,'[SOURCE EQUIVALENT]'))
   tline = fgetl(fid);
   
   if ~isempty(strfind(tline,'FORMAT_SOURCEEQUIVALENT=NodeID,Voltage,PhaseAngle,PositiveSequenceResistance,PositiveSequenceReactance,ZeroSequenceResistance,ZeroSequenceReactance,Configuration,KVLLdesired,ImpedanceUnit'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[LOAD EQUIVALENT]')) 
             ySOURCEEQUIVALENT=tline;
             fprintf(fid2,'%s\n',ySOURCEEQUIVALENT);
             tline = fgetl(fid);
       end
   end    
% HydroOttawa system
%    if ~isempty(strfind(tline,'FORMAT_SOURCEEQUIVALENT=NodeID,Voltage,OperatingAngle1,OperatingAngle2,OperatingAngle3,PositiveSequenceResistance,PositiveSequenceReactance,ZeroSequenceResistance,ZeroSequenceReactance,Configuration,OperatingVoltage1,OperatingVoltage2,OperatingVoltage3,ImpedanceUnit'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[LOAD EQUIVALENT]')) 
%              ySOURCEEQUIVALENT=tline;
%              fprintf(fid2,'%s\n\n',ySOURCEEQUIVALENT);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end
% 
% LOAD EQUIVALENT data
fid2 = fopen([FileName '_LOADEQUIVALENT.txt'], 'w');
if ~isempty(strfind(tline,'[LOAD EQUIVALENT]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_LOADEQUIVALENT=NodeID,Format,Value1A,Value1B,Value1C,Value2A,Value2B,Value2C'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[LINE CONFIGURATION]')) 
             yLOADEQUIVALENT=tline;
             fprintf(fid2,'%s\r\n',yLOADEQUIVALENT);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% Line Configuration data
fid2 = fopen([FileName '_LINECONFIGURATION.txt'], 'w');
if ~isempty(strfind(tline,'[LINE CONFIGURATION]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_LINECONFIGURATION=SectionID,DeviceNumber,LineCableID,Length,Overhead,NumberOfCableInParallel,ConnectionStatus,CoordX,CoordY'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SECTION]')) 
             yLINECONFIGURATION=tline;
             fprintf(fid2,'%s\r\n',yLINECONFIGURATION);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% fid2 = fopen([FileName '_OVERHEADLINESETTING.txt'], 'w');
% if ~isempty(strfind(tline,'[OVERHEADLINE SETTING]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_OVERHEADLINESETTING=SectionID,DeviceNumber,LineCableID,Length,ConnectionStatus,CoordX,CoordY,HarmonicModel,FlowConstraintActive,FlowConstraintUnit,MaximumFlow,SeriesCompensationActive,MaxReactanceMultiplier,SeriesCompensationCost'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[UNDERGROUNDLINE SETTING]')) 
%              yOVERHEADLINESETTING=tline;
%              fprintf(fid2,'%s\r\n',yOVERHEADLINESETTING);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% fid2 = fopen([FileName '_UNDERGROUNDLINESETTING.txt'], 'w');
% if ~isempty(strfind(tline,'[UNDERGROUNDLINE SETTING]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_UNDERGROUNDLINESETTING=SectionID,DeviceNumber,LineCableID,Length,NumberOfCableInParallel,Amps,Amps_1,Amps_2,Amps_3,Amps_4,ConnectionStatus,CoordX,CoordY,HarmonicModel,FlowConstraintActive,FlowConstraintUnit,MaximumFlow'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SECTION]')) 
%              yUNDERGROUNDLINESETTING=tline;
%              fprintf(fid2,'%s\r\n',yUNDERGROUNDLINESETTING);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end

% SECTION data
fid2 = fopen([FileName '_SECTION.txt'], 'w');
if ~isempty(strfind(tline,'[SECTION]'))
   tline = fgetl(fid);
    if ~isempty(strfind(tline,'FORMAT_SECTION=SectionID,FromNodeID,ToNodeID,Phase,ZoneID,SubNetworkId,EnvironmentID'))
        tline = fgetl(fid);
        if ~isempty(strfind(tline,'FORMAT_FEEDER=NetworkID,HeadNodeID,CoordSet,Year,Description,Color,LoadFactor,Group1,Group2,TagText,TagProperties,TagDeltaX,TagDeltaY,TagAngle,TagAlignment,TagBorder,TagBackground,TagTextColor,TagBorderColor,TagBackgroundColor,TagLocation,TagFont,TagTextSize,Version,EnvironmentID'))
            tline = fgetl(fid);
%             if ~isempty(strfind(tline,'FEEDER=22M32S,22M32S-11649310,1,0,,0,1.000000,_HV SOURCE_REG,230KV,,NULL,,,,,,,,,,,,,,,-1,0'))
%                 tline = fgetl(fid);
                while isempty(strfind(tline,'[SWITCH SETTING]')) 
                    ySECTION=tline;
                    fprintf(fid2,'%s\r\n',ySECTION);
                    tline = fgetl(fid);
                end
%             end
        end 

%            if ~isempty(strfind(tline,'FORMAT_FEEDER=NetworkID,HeadNodeID,CoordSet,Year,Description,Color,LoadFactor,Group1,Group2,Group3,TagText,TagProperties,TagDeltaX,TagDeltaY,TagAngle,TagAlignment,TagBorder,TagBackground,TagTextColor,TagBorderColor,TagBackgroundColor,TagLocation,TagFont,TagTextSize,TagOffset,Version,EnvironmentID'))
%               tline = fgetl(fid);
%              if ~isempty(strfind(tline,'FEEDER=22M32S,22M32S-11649310,1,0,,0,1.000000,_HV SOURCE_REG,230KV,,NULL,,,,,,,,,,,,,,,-1,0'))
%                   tline = fgetl(fid);
%                    while isempty(strfind(tline,'[SUBNETWORKS]')) 
%                          ySECTION=tline;
%                          fprintf(fid2,'%s\r\n',ySECTION);
%                          tline = fgetl(fid);
%                    end
%              end
%            end 
   end    
   fclose(fid2);
end

% % SUBNETWORKS data
% fid2 = fopen([FileName '_SUBNETWORKS.txt'], 'w');
% if ~isempty(strfind(tline,'[SUBNETWORKS]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_SUBNETWORKS=SubNetID,Angle,X,Y,Height,Length,ParentSubNetID,SymbolID,TagText,TagProperties,TagDeltaX,TagDeltaY,TagAngle,TagAlignment,TagBorder,TagBackground,TagTextColor,TagBorderColor,TagBackgroundColor,TagLocation,TagFont,TagTextSize,TagOffset,SubNetTypeId,Version,CoordSet'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SUBNETWORK CONNECTIONS]')) 
%              ySUBNETWORKS=tline;
%              fprintf(fid2,'%s\r\n',ySUBNETWORKS);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% % SUBNETWORK CONNECTIONS data
% fid2 = fopen([FileName '_SUBNETWORKCONNECTIONS.txt'], 'w');
% if ~isempty(strfind(tline,'[SUBNETWORK CONNECTIONS]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_SUBNETWORKCONNECTIONS=SubNetID,NodeID,ConnectorCoordX,ConnectorCoordY,ConnectorIndex,SymbolConnectorIndex,Description'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[TRANSFORMER SETTING]')) 
%              ySUBNETWORKCONNECTIONS=tline;
%              fprintf(fid2,'%s\r\n',ySUBNETWORKCONNECTIONS);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% 
% % TRANSFORMERSETTING data
% fid2 = fopen([FileName '_TRANSFORMERSETTING.txt'], 'w');
% if ~isempty(strfind(tline,'[TRANSFORMER SETTING]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_TRANSFORMERSETTING=SectionID,Location,EqID,DeviceNumber,CoordX,CoordY,Conn,PrimTap,SecondaryTap,RgPrim,XgPrim,RgSec,XgSec,ODPrimPh,PrimaryBaseVoltage,SecondaryBaseVoltage,FromNodeID,SettingOption,SetPoint,ControlType,LowerBandWidth,UpperBandWidth,TapLocation,InitialTapPosition,InitialTapPositionMode,Tap,MaxBuck,MaxBoost,CT,PT,Rset,Xset,FirstHouseHigh,FirstHouseLow,PhaseON,AtSectionID,MasterID,FaultIndicator,DiversityFactorA,DiversityFactorB,DiversityFactorC,PhaseShiftType,GammaPhaseShift,DemandType,Value1A,Value2A,Value1B,Value2B,Value1C,Value2C,Value1ABC,Value2ABC,ConnectionStatus,Reversible'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SWITCH SETTING]')) 
%              yTRANSFORMERSETTING=tline;
%              fprintf(fid2,'%s\r\n',yTRANSFORMERSETTING);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end


% SWITCH SETTING data
fid2 = fopen([FileName '_SWITCHSETTING.txt'], 'w');
if ~isempty(strfind(tline,'[SWITCH SETTING]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'SectionID,Location,EqID,DeviceNumber,EqPhase,CoordX,CoordY,ClosedPhase,Locked,RC,NStatus,DemandType,Value1A,Value2A,Value1B,Value2B,Value1C,Value2C,Value1ABC,Value2ABC,TCCID,Desc,PhPickup,GrdPickup,Alternate,PhAltPickup,GrdAltPickup,FromNodeID,FaultIndicator,Automated,SensorMode,Strategic,RestorationMode,ConnectionStatus,ByPassOnRestoration'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[BREAKER SETTING]')) 
             ySWITCHSETTING=tline;
             fprintf(fid2,'%s\r\n',ySWITCHSETTING);
             tline = fgetl(fid);
       end
   end  
%    if ~isempty(strfind(tline,'FORMAT_SWITCHSETTING=SectionID,Location,EqID,DeviceNumber,CoordX,CoordY,ClosedPhase,Locked,RC,NStatus,DemandType,Value1A,Value2A,Value1B,Value2B,Value1C,Value2C,Value1ABC,Value2ABC,PhPickup,GrdPickup,Alternate,PhAltPickup,GrdAltPickup,FromNodeID,FaultIndicator,Automated,SensorMode,Strategic,RestorationMode,ConnectionStatus,ByPassOnRestoration,Reversible'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[BREAKER SETTING]')) 
%              ySWITCHSETTING=tline;
%              fprintf(fid2,'%s\r\n',ySWITCHSETTING);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% BREAKER SETTING data
fid2 = fopen([FileName '_BREAKERSETTING.txt'], 'w');
if ~isempty(strfind(tline,'[BREAKER SETTING]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_BREAKERSETTING=SectionID,Location,EqID,DeviceNumber,EqPhase,CoordX,CoordY,ClosedPhase,Locked,RC,NStatus,DemandType,Value1A,Value2A,Value1B,Value2B,Value1C,Value2C,Value1ABC,Value2ABC,TCCID,Desc,PhPickup,GrdPickup,Alternate,PhAltPickup,GrdAltPickup,FromNodeID,EnableReclosing,FaultIndicator,EnableFuseSaving,MinRatedCurrentForFuseSaving,Automated,SensorMode,Strategic,RestorationMode,ConnectionStatus,ByPassOnRestoration'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SHUNT CAPACITOR SETTING]')) 
             yBREAKERSETTING=tline;
             fprintf(fid2,'%s\r\n',yBREAKERSETTING);
             tline = fgetl(fid);
       end
   end  
%    if ~isempty(strfind(tline,'FORMAT_BREAKERSETTING=SectionID,Location,EqID,DeviceNumber,CoordX,CoordY,ClosedPhase,Locked,RC,NStatus,DemandType,Value1A,Value2A,Value1B,Value2B,Value1C,Value2C,Value1ABC,Value2ABC,TCCID,PhPickup,GrdPickup,Alternate,PhAltPickup,GrdAltPickup,FromNodeID,EnableReclosing,FaultIndicator,EnableFuseSaving,MinRatedCurrentForFuseSaving,Automated,SensorMode,Strategic,RestorationMode,ConnectionStatus,ByPassOnRestoration,Speed,Reversible'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[FUSE SETTING]')) 
%              yBREAKERSETTING=tline;
%              fprintf(fid2,'%s\r\n',yBREAKERSETTING);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% SHUNT CAPACITOR SETTING data
fid2 = fopen([FileName '_SHUNTCAPACITORSETTING.txt'], 'w');
if ~isempty(strfind(tline,'[SHUNT CAPACITOR SETTING]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SHUNTCAPACITORSETTING=SectionID,DeviceNumber,Location,Connection,Phase,KVAR,ThreePhaseKVAR,Losses,ThreePhaseLosses,KVLN,Control,ONSetting,OFFSetting,ShuntCapacitorID,ControllingPhase,ConnectionStatus,SwitchedCapacitorStatus'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[INTERMEDIATE NODES]')) 
             ySHUNTCAPACITORSETTING=tline;
             fprintf(fid2,'%s\r\n',ySHUNTCAPACITORSETTING);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end


% INTERMEDIATE NODES data; 
% no 'GetCYMData' for this one, yet; 
% does not seem to be necessary
fid2 = fopen([FileName '_INTERMEDIATENODES.txt'], 'w');
if ~isempty(strfind(tline,'[INTERMEDIATE NODES]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_INTERMEDIATENODE=SectionID,SeqNumber,CoordX,CoordY,IsBreakPoint,BreakPointLocation'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'    '))  
             yINTERMEDIATENODES=tline;
             fprintf(fid2,'%s\r\n',yINTERMEDIATENODES);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% 
% % RECLOSER SETTING data
% fid2 = fopen([FileName '_RECLOSERSETTING.txt'], 'w');
% if ~isempty(strfind(tline,'[RECLOSER SETTING]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_RECLOSERSETTING=SectionID,Location,EqID,DeviceNumber,CoordX,CoordY,ClosedPhase,Locked,RC,NStatus,DemandType,Value1A,Value2A,Value1B,Value2B,Value1C,Value2C,Value1ABC,Value2ABC,TCCID,PhPickup,GrdPickup,Alternate,PhAltPickup,GrdAltPickup,FromNodeID,EnableReclosing,FaultIndicator,EnableFuseSaving,MinRatedCurrentForFuseSaving,Automated,SensorMode,Strategic,RestorationMode,ConnectionStatus,ByPassOnRestoration,Reversible,TCCSettingsSelection'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[INTERMEDIATE NODES]')) 
%              yRECLOSERSETTING=tline;
%              fprintf(fid2,'%s\r\n',yRECLOSERSETTING);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% % INTERMEDIATE NODES data
% fid2 = fopen([FileName '_INTERMEDIATENODES.txt'], 'w');
% if ~isempty(strfind(tline,'[INTERMEDIATE NODES]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_INTERMEDIATENODE=SectionID,SeqNumber,CoordX,CoordY,IsBreakPoint,BreakPointLocation'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[DEVICETAG]')) 
%              yINTERMEDIATENODES=tline;
%              fprintf(fid2,'%s\r\n',yINTERMEDIATENODES);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% % DEVICETAG data
% fid2 = fopen([FileName '_DEVICETAG.txt'], 'w');
% if ~isempty(strfind(tline,'[DEVICETAG]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_DEVICETAG=DeviceNumber,DeviceType,TagText,TagProperties,TagDeltaX,TagDeltaY,TagAngle,TagAlignment,TagBorder,TagBackground,TagTextColor,TagBorderColor,TagBackgroundColor,TagLocation,TagFont,TagTextSize,TagOffset'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'AUTO TAP CHANGING EXTENSION')) 
%              yDEVICETAG=tline;
%              fprintf(fid2,'%s\r\n',yDEVICETAG);
%              tline = fgetl(fid);
%              if ~ischar(tline),   break
%              end
%        end
%    end    
%    fclose(fid2);
% end
% 
% % DEVICETAG data
% fid2 = fopen([FileName '_AUTOTAPCHANGINGEXTENSION.txt'], 'w');
% if ~isempty(strfind(tline,'[AUTO TAP CHANGING EXTENSION]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_AUTOTAPCHANGINGEXTST=DeviceNumber,DeviceType,LTCIndex,TapChangingAlgo,ActivationDelayFirstTap,ActivationDelaySubTap,MechanismDelay,ControlOnTapSide,WantIntegration,UseInverseTimeDelayChar'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'    ')) 
%              yAUTOTAPCHANGINGEXTENSION=tline;
%              fprintf(fid2,'%s\r\n',yAUTOTAPCHANGINGEXTENSION);
%              tline = fgetl(fid);
%              if ~ischar(tline),   break
%              end
%        end
%    end    
%    fclose(fid2);
% end

ExecutionFlag = 1;

fclose(fid);



    