function[ExecutionFlag] = GetCYMData_equipment_SCE(FilePath)


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


yCONCENTRICNEUTRALCABLE=[];
%ySHIELDEDCABLE=[];
%yUNSHIELDEDCABLE=[];
yCONDUCTOR=[];
ySPACINGTABLEFORLINE=[];
yDOUBLECIRCUITSPACING=[];
yLINE=[];
yLINEUNBALANCED=[];
yINDUCTIONMOTOR=[];
ySYNCHRONOUSMOTOR=[];
yINDUCTIONGENERATOR=[];
ySYNCHRONOUSGENERATOR=[];
yELECTRONICCONVERTERGENERATOR=[];
yREGULATOR=[];
ySWITCH=[];
yBREAKER=[];
yFUSE=[];
yLVCB=[];
yMISCELLANEOUS=[];
yRECLOSER=[];
ySECTIONALIZER=[];
ySUBSTATION=[];
yTRANSFORMER=[];
yTHREEWINDINGTRANSFORMER=[];
yAUTOTRANSFORMER=[];
yTHREEWINDINGAUTOTRANSFORMER=[];
ySERIECAPACITOR=[];
ySERIEREACTOR=[];
ySHUNTCAPACITOR=[];
ySHUNTREACTOR=[];
yARCFURNACE=[];
yCTYPEFILTER=[];
yDOUBLETUNEDFILTER=[];
yHIGHPASSFILTER=[];
yIDEALCONVERTER=[];
yNONIDEALCONVERTER=[];
yFREQUENCYSOURCE=[];
ySINGLETUNEDFILTER=[];
yGROUNDINGTRANSFORMER=[];
yWECS=[];
yINDUCTIONMACHINEEQCIRCUIT=[];
ySYNCHRONOUSMACHINEEQCIRCUIT=[];
ySYNCHRONOUSMACHINEEXTSTAB=[];
    ySYNCHRONOUSMACHINEEXTHARMO=[];
yWINDMODEL=[];
%yWINDMODELPOINTS=[];
%yDEFAULTEQUIPMENT=[];
yINSOLATIONMODEL=[];
%yINSOLATIONMODELPOINTBLOCK=[];
yPHOTOVOLTAIC=[];
ySOFC=[];
yMICROTURBINE=[];
yNETWORKPROTECTOR=[];
ySVC=[];
yGENERATORCOSTCURVEMODEL=[];
yBUSWAY=[];
yPHASESHIFTERTRANSFORMER=[];
yGENERATIONCURVEMODEL=[];
yLOADCURVEMODEL=[];
yMOTORCURVEMODEL=[];
yVARIABLEFREQUENCYDRIVE=[];



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

% CONCENTRIC NEUTRAL CABLE data
fid2 = fopen([FileName '_CONCENTRICNEUTRALCABLE.txt'], 'w');
tline = fgetl(fid);
while isempty(strfind(tline,'[SHIELDED CABLE]'))
      yCONCENTRICNEUTRALCABLE=tline;
      fprintf(fid2,'%s\r\n',yCONCENTRICNEUTRALCABLE);
      tline = fgetl(fid);
end
fclose(fid2);

% SHIELDED CABLE data
fid2 = fopen([FileName '_SHIELDEDCABLE.txt'], 'w');
tline = fgetl(fid);
while isempty(strfind(tline,'[CONDUCTOR]'))
      ySHIELDEDCABLE=tline;
      fprintf(fid2,'%s\r\n',ySHIELDEDCABLE);
      tline = fgetl(fid);
end
fclose(fid2);


% CONDUCTOR data
fid2 = fopen([FileName '_CONDUCTOR.txt'], 'w');
if ~isempty(strfind(tline,'[CONDUCTOR]'))
        disp '4'
   tline = fgetl(fid);
   
   if ~isempty(strfind(tline,'FORMAT_CONDUCTOR=ID,Diameter,GMR,R25,R50,Amps,Amps_2,WithstandRating,FailRate,TmpFailRate,OutageTime,Material,CodeWord,Size,Standard,ConstructionType,InsideDiameter,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SPACING TABLE FOR LINE]')) 
         yCONDUCTOR=tline;
         fprintf(fid2,'%s\r\n',yCONDUCTOR);
         tline = fgetl(fid);
       end
   end    
%    if ~isempty(strfind(tline,'FORMAT_CONDUCTOR=ID,Diameter,GMR,R25,R50,Amps,Amps_1,Amps_2,Amps_3,Amps_4,WithstandRating,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,Material,CodeWord,ConstructionType,InsideDiameter,FirstResistanceDC,SecondResistanceDC,Size_mm2,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SPACING TABLE FOR LINE]')) 
%          yCONDUCTOR=tline;
%          fprintf(fid2,'%s\r\n',yCONDUCTOR);
%          tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end
 

% SPACING TABLE FOR LINE data
fid2 = fopen([FileName '_SPACINGTABLEFORLINE.txt'], 'w');
if ~isempty(strfind(tline,'[SPACING TABLE FOR LINE]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SPACINGTABLEFORLINE=ID,GMDPh-Ph,GMDPh-N,AvgPhCondHeight,AvgNeutralHeight,PosOfCond1_X,PosOfCond1_Y,PosOfCond2_X,PosOfCond2_Y,PosOfCond3_X,PosOfCond3_Y,PosOfNeutralCond_X,PosOfNeutralCond_Y,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[LINE]')) 
         ySPACINGTABLEFORLINE=tline;
         fprintf(fid2,'%s\r\n',ySPACINGTABLEFORLINE);
         tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end
 
% % SPACING TABLE FOR LINE data
% fid2 = fopen([FileName '_DOUBLECIRCUITSPACING.txt'], 'w');
% if ~isempty(strfind(tline,'[DOUBLECIRCUITSPACING]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_DOUBLECIRCUITSPACING=ID,PosOfCond1_X,PosOfCond1_Y,PosOfCond2_X,PosOfCond2_Y,PosOfCond3_X,PosOfCond3_Y,PosOfNeutralCond_X,PosOfNeutralCond_Y,PosOfCond1_C2_X,PosOfCond1_C2_Y,PosOfCond2_C2_X,PosOfCond2_C2_Y,PosOfCond3_C2_X,PosOfCond3_C2_Y,PosOfNeutralCond_N2_X,PosOfNeutralCond_N2_Y,BundleDistance,NBPhasesPerCircuit,NBConductorsPerPhase,NBNeutrals,TowerType,DistanceA,DistanceB,DistanceC,DistanceD,DistanceE,ConductorStatusN1,ConductorStatusN2,FootingResistanceN1,FootingResistanceN2,TowerSpanN1,TowerSpanN2,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[LINE]')) 
%          yDOUBLECIRCUITSPACING=tline;
%          fprintf(fid2,'%s\r\n',yDOUBLECIRCUITSPACING);
%          tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end

% LINE data
fid2 = fopen([FileName '_LINE.txt'], 'w');
if ~isempty(strfind(tline,'[LINE]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_LINE=ID,PhaseCondID,NeutralCondID,SpacingID,R1,R0,X1,X0,B1,B0,Amps,Amps_2,LockImpedance,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[LINE UNBALANCED]')) 
             yLINE=tline;
             fprintf(fid2,'%s\r\n',yLINE);
             tline = fgetl(fid);
       end
   end 
%    if ~isempty(strfind(tline,'FORMAT_LINE=ID,PhaseCondID,NeutralCondID,SpacingID,R1,R0,X1,X0,B1,B0,Amps,Amps_1,Amps_2,Amps_3,Amps_4,LockImpedance,PositiveSequenceShuntConductance,ZeroSequenceShuntConductance,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[LINE UNBALANCED]')) 
%              yLINE=tline;
%              fprintf(fid2,'%s\r\n',yLINE);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% LINE UNBALANCED data
fid2 = fopen([FileName '_LINEUNBALANCED.txt'], 'w');
if ~isempty(strfind(tline,'[LINE UNBALANCED]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_LINEUNBALANCED=ID,CondID_A,CondID_B,CondID_C,CondID_N,SpacingID,Ra,Rb,Rc,Rm,Xa,Xb,Xc,Xm,Ba,Bb,Bc,Bm,AmpsA,AmpsB,AmpsC,AmpsA_2,AmpsB_2,AmpsC_2,LockImpedance,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[INDUCTION MOTOR]')) 
             yLINEUNBALANCED=tline;
             fprintf(fid2,'%s\r\n',yLINEUNBALANCED);
             tline = fgetl(fid);
       end
   end    

%    if ~isempty(strfind(tline,'FORMAT_LINEUNBALANCED=ID,CondID_A,CondID_B,CondID_C,CondID_N,SpacingID,Ra,Rb,Rc,Xa,Xb,Xc,Ba,Bb,Bc,AmpsA,AmpsB,AmpsC,AmpsA_1,AmpsB_1,AmpsC_1,AmpsA_2,AmpsB_2,AmpsC_2,AmpsA_3,AmpsB_3,AmpsC_3,AmpsA_4,AmpsB_4,AmpsC_4,LockImpedance,ShuntConductanceA,ShuntConductanceB,ShuntConductanceC,MutualResistanceAB,MutualResistanceBC,MutualResistanceCA,MutualReactanceAB,MutualReactanceBC,MutualReactanceCA,MutualShuntSusceptanceAB,MutualShuntSusceptanceBC,MutualShuntSusceptanceCA,MutualShuntConductanceAB,MutualShuntConductanceBC,MutualShuntConductanceCA,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[INDUCTION MOTOR]')) 
%              yLINEUNBALANCED=tline;
%              fprintf(fid2,'%s\r\n',yLINEUNBALANCED);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% INDUCTION MOTOR data
fid2 = fopen([FileName '_INDUCTIONMOTOR.txt'], 'w');
if ~isempty(strfind(tline,'[INDUCTION MOTOR]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_INDUCTIONMOTOR=ID,RatedPower,RatedVoltageKVLL,Efficiency,FullLoadPF,Type,RatedSpeed,SubTransientResistance,SubTransientReactance,ANSIMotorGroup,LockedRotorPF,NEMA,KVA_HPRatio,SymbolID,ComputeMode,ImpedanceUnit,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SYNCHRONOUS MOTOR]')) 
             yINDUCTIONMOTOR=tline;
             fprintf(fid2,'%s\r\n',yINDUCTIONMOTOR);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% SYNCHRONOUS MOTOR data
fid2 = fopen([FileName '_SYNCHRONOUSMOTOR.txt'], 'w');
if ~isempty(strfind(tline,'SYNCHRONOUS MOTOR'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SYNCHRONOUSMOTOR=ID,RatedVoltageKVLL,RatedPower,Efficiency,FullLoadPF,SynchronousSpeed,NumberOfPoles,Connection,SubTransientResistance,SubTransientReactance,XdSaturated,ZeroSequenceResistance,ZeroSequenceReactance,GroundingResistance,GroundingReactance,SymbolID,ImpedanceUnit,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[INDUCTION GENERATOR]')) 
             ySYNCHRONOUSMOTOR=tline;
             fprintf(fid2,'%s\r\n',ySYNCHRONOUSMOTOR);
             tline = fgetl(fid);
       end
   end 
%    if ~isempty(strfind(tline,'FORMAT_SYNCHRONOUSMOTOR=ID,RatedVoltageKVLL,RatedPower,Efficiency,FullLoadPF,SynchronousSpeed,NumberOfPoles,Connection,SubTransientResistance,SubTransientReactance,XdSaturated,ZeroSequenceResistance,ZeroSequenceReactance,GroundingResistance,GroundingReactance,SymbolID,ImpedanceUnit,NegativeSequenceResistance,NegativeSequenceReactance,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[INDUCTION GENERATOR]')) 
%              ySYNCHRONOUSMOTOR=tline;
%              fprintf(fid2,'%s\r\n',ySYNCHRONOUSMOTOR);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end
% INDUCTION GENERATOR data
fid2 = fopen([FileName '_INDUCTIONGENERATOR.txt'], 'w');
if ~isempty(strfind(tline,'[INDUCTION GENERATOR]'))
   tline = fgetl(fid);
    if ~isempty(strfind(tline,'FORMAT_INDUCTIONGENERATOR=ID,KVA,KVLL,ActiveGeneration,PF,RatedSpeed,ANSIMotorGroup,SubTransientResistance,SubTransientReactance,SymbolID,AutoComputeFromEqCircuit,ImpedanceUnit,Efficiency,Comments'))
       tline = fgetl(fid);
           while isempty(strfind(tline,'[SYNCHRONOUS GENERATOR]')) 
                 yINDUCTIONGENERATOR=tline;
                 fprintf(fid2,'%s\r\n',yINDUCTIONGENERATOR);
                 tline = fgetl(fid);
           end
   end    
   fclose(fid2);
end
% SYNCHRONOUS GENERATOR data
fid2 = fopen([FileName '_SYNCHRONOUSGENERATOR.txt'], 'w');
if ~isempty(strfind(tline,'[SYNCHRONOUS GENERATOR]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SYNCHRONOUSGENERATOR=ID,KVA,KVLL,ActiveGeneration,PF,MaxKVAR,MinKVAR,ConnectionConfiguration,R1,X1,R0,X0,Rg,Xg,TransientResistance,TransientReactance,SubtransientResistance,SubtransientReactance,NumberOfPoles,SymbolID,ImpedanceUnit,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[ELECTRONIC CONVERTER GENERATOR]')) 
             ySYNCHRONOUSGENERATOR=tline;
             fprintf(fid2,'%s\r\n',ySYNCHRONOUSGENERATOR);
             tline = fgetl(fid);
       end
   end    
%    if ~isempty(strfind(tline,'FORMAT_SYNCHRONOUSGENERATOR=ID,KVA,KVLL,ActiveGeneration,PF,MaxKVAR,MinKVAR,ConnectionConfiguration,R1,X1,R0,X0,Rg,Xg,TransientResistance,TransientReactance,SubtransientResistance,SubtransientReactance,NumberOfPoles,SymbolID,ImpedanceUnit,FixedQLimits,NegativeSequenceResistance,NegativeSequenceReactance,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[ELECTRONIC CONVERTER GENERATOR]')) 
%              ySYNCHRONOUSGENERATOR=tline;
%              fprintf(fid2,'%s\r\n',ySYNCHRONOUSGENERATOR);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% ELECTRONIC CONVERTER GENERATOR data
fid2 = fopen([FileName '_ELECTRONICCONVERTERGENERATOR.txt'], 'w');
if ~isempty(strfind(tline,'[ELECTRONIC CONVERTER GENERATOR]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_ELECTRONICCONVERTERGENERATOR=ID,KVA,KVLL,ActiveGeneration,PF,FaultContribution,ANSIMotorGroup,Converter,SymbolID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[REGULATOR]')) 
             yELECTRONICCONVERTERGENERATOR=tline;
             fprintf(fid2,'%s\r\n',yELECTRONICCONVERTERGENERATOR);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% REGULATOR data
fid2 = fopen([FileName '_REGULATOR.txt'], 'w');
if ~isempty(strfind(tline,'[REGULATOR]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_REGULATOR=ID,Type,KVA,KVA_2,KVLN,MaxBuck,MaxBoost,Taps,Bandwidth,CT,PT,Reversible,SymbolID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SWITCH]')) 
             yREGULATOR=tline;
             fprintf(fid2,'%s\r\n',yREGULATOR);
             tline = fgetl(fid);
       end
   end  
%    if ~isempty(strfind(tline,'FORMAT_REGULATOR=ID,Type,KVA,KVA_1,KVA_2,KVA_3,KVA_4,KVLN,MaxBuck,MaxBoost,Taps,Bandwidth,CT,PT,Reversible,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,ConstructionType,EnableBonusRating,BonusRating1,BonusRating2,BonusRating3,BonusRating4,RegulationRange1,RegulationRange2,RegulationRange3,RegulationRange4,MaxRatingCurrent,SymbolID,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SWITCH]')) 
%              yREGULATOR=tline;
%              fprintf(fid2,'%s\r\n',yREGULATOR);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% SWITCH data
fid2 = fopen([FileName '_SWITCH.txt'], 'w');
if ~isempty(strfind(tline,'[SWITCH]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SWITCH=ID,Amps,Amps_2,KVLL,Reversible,FailRate,TmpFailRate,OutageTime,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,RemoteControlled,Automated,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[BREAKER]')) 
             ySWITCH=tline;
             fprintf(fid2,'%s\r\n',ySWITCH);
             tline = fgetl(fid);
       end
   end    
%    if ~isempty(strfind(tline,'FORMAT_SWITCH=ID,Amps,Amps_1,Amps_2,Amps_3,Amps_4,KVLL,Reversible,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,RemoteControlled,Automated,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[BREAKER]')) 
%              ySWITCH=tline;
%              fprintf(fid2,'%s\r\n',ySWITCH);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% BREAKER data
fid2 = fopen([FileName '_BREAKER.txt'], 'w');
if ~isempty(strfind(tline,'[BREAKER]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_BREAKER=ID,Amps,Amps_2,KVLL,Reversible,InterruptingRating,FailRate,TmpFailRate,OutageTime,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,SinglePhaseTripping,RemoteControlled,Automated,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[FUSE]')) 
             yBREAKER=tline;
             fprintf(fid2,'%s\r\n',yBREAKER);
             tline = fgetl(fid);
       end
   end 
%    if ~isempty(strfind(tline,'FORMAT_BREAKER=ID,Amps,Amps_1,Amps_2,Amps_3,Amps_4,KVLL,Reversible,InterruptingRating,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,SinglePhaseTripping,RemoteControlled,Automated,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[FUSE]')) 
%              yBREAKER=tline;
%              fprintf(fid2,'%s\r\n',yBREAKER);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% FUSE data
fid2 = fopen([FileName '_FUSE.txt'], 'w');
if ~isempty(strfind(tline,'[FUSE]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_FUSE=ID,Amps,Amps_2,KVLL,Reversible,InterruptingRating,FailRate,TmpFailRate,OutageTime,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[LVCB]')) 
             yFUSE=tline;
             fprintf(fid2,'%s\r\n',yFUSE);
             tline = fgetl(fid);
       end
   end 
%    if ~isempty(strfind(tline,'FORMAT_FUSE=ID,Amps,Amps_1,Amps_2,Amps_3,Amps_4,KVLL,Reversible,InterruptingRating,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,Comments,Manufacturer,Model,TCCRating'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[LVCB]')) 
%              yFUSE=tline;
%              fprintf(fid2,'%s\r\n',yFUSE);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% LVCB data
fid2 = fopen([FileName '_LVCB.txt'], 'w');
if ~isempty(strfind(tline,'[LVCB]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_LVCB=ID,Amps,Amps_2,KVLL,Reversible,InterruptingRating,FailRate,TmpFailRate,OutageTime,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,SinglePhaseTripping,RemoteControlled,Automated,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[MISCELLANEOUS]')) 
             yLVCB=tline;
             fprintf(fid2,'%s\r\n',yLVCB);
             tline = fgetl(fid);
       end
   end 
%    if ~isempty(strfind(tline,'FORMAT_LVCB=ID,Amps,Amps_1,Amps_2,Amps_3,Amps_4,KVLL,Reversible,InterruptingRating,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,SinglePhaseTripping,RemoteControlled,Automated,Comments,LVCBType,Manufacturer,Model'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[MISCELLANEOUS]')) 
%              yLVCB=tline;
%              fprintf(fid2,'%s\r\n',yLVCB);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end
% MISCELLANEOUS data
fid2 = fopen([FileName '_MISCELLANEOUS.txt'], 'w');
if ~isempty(strfind(tline,'[MISCELLANEOUS]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_MISCELLANEOUS=ID,Amps,Amps_2,KVLL,SymbolOpenID,SymbolCloseID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[RECLOSER]')) 
             yMISCELLANEOUS=tline;
             fprintf(fid2,'%s\r\n',yMISCELLANEOUS);
             tline = fgetl(fid);
       end
   end  
%    if ~isempty(strfind(tline,'FORMAT_MISCELLANEOUS=ID,Amps,Amps_1,Amps_2,Amps_3,Amps_4,KVLL,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,SymbolOpenID,SymbolCloseID,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[RECLOSER]')) 
%              yMISCELLANEOUS=tline;
%              fprintf(fid2,'%s\r\n',yMISCELLANEOUS);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% RECLOSER data
fid2 = fopen([FileName '_RECLOSER.txt'], 'w');
if ~isempty(strfind(tline,'[RECLOSER]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_RECLOSER=ID,Amps,Amps_2,KVLL,Reversible,InterruptingRating,FailRate,TmpFailRate,OutageTime,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,SinglePhaseTripping,RemoteControlled,Automated,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[RELAY]')) 
         yRECLOSER=tline;
         fprintf(fid2,'%s\r\n',yRECLOSER);
         tline = fgetl(fid);
       end
   end    
%    if ~isempty(strfind(tline,'FORMAT_RECLOSER=ID,Amps,Amps_1,Amps_2,Amps_3,Amps_4,KVLL,Reversible,InterruptingRating,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,SinglePhaseTripping,RemoteControlled,Automated,Comments,RecloserType,ControlType,Model'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SECTIONALIZER]')) 
%          yRECLOSER=tline;
%          fprintf(fid2,'%s\r\n',yRECLOSER);
%          tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% RELAY data
fid2 = fopen([FileName '_RELAY.txt'], 'w');
if ~isempty(strfind(tline,'[RELAY]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_RELAY=ID,Amps,Amps_2,KVLL,Reversible,InterruptingRating,FailRate,TmpFailRate,OutageTime,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,SinglePhaseTripping,RemoteControlled,Automated,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SECTIONALIZER]')) 
             yRELAY=tline;
             fprintf(fid2,'%s\r\n',yRELAY);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% SECTIONALIZER data
fid2 = fopen([FileName '_SECTIONALIZER.txt'], 'w');
if ~isempty(strfind(tline,'[SECTIONALIZER]'))
   tline = fgetl(fid);
    if ~isempty(strfind(tline,'FORMAT_SECTIONALIZER=ID,Amps,Amps_2,KVLL,Reversible,FailRate,TmpFailRate,OutageTime,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,RemoteControlled,Automated,InterruptingRating,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SUBSTATION]')) 
             ySECTIONALIZER=tline;
             fprintf(fid2,'%s\r\n',ySECTIONALIZER);
             tline = fgetl(fid);
       end
   end    

%    if ~isempty(strfind(tline,'FORMAT_SECTIONALIZER=ID,Amps,Amps_1,Amps_2,Amps_3,Amps_4,KVLL,Reversible,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,RemoteControlled,Automated,InterruptingRating,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SUBSTATION]')) 
%              ySECTIONALIZER=tline;
%              fprintf(fid2,'%s\r\n',ySECTIONALIZER);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% SUBSTATION data
fid2 = fopen([FileName '_SUBSTATION.txt'], 'w');
if ~isempty(strfind(tline,'[SUBSTATION]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SUBSTATION=ID,MVA,MVA_2,KVLL,KVLLdesired,R1,X1,R0,X0,Conn,PhaseAngle,PrimaryEquivalentType,SubEqVal1,SubEqVal2,SubEqVal3,SubEqVal4,SubPrimayLLVoltage,SecondaryFaultReactance,TxfoConnection,HarmonicEnveloppe,ImpedanceUnit,BranchID_1,PrimProtDevID_1,PrimProtDevNum_1,TransformerID_1,TransformerNum_1,SubXs_1,SecProtDevID_1,SecProtDevNum_1,BranchStatus_1,BranchID_2,PrimProtDevID_2,PrimProtDevNum_2,TransformerID_2,TransformerNum_2,SubXs_2,SecProtDevID_2,SecProtDevNum_2,BranchStatus_2,BranchID_3,PrimProtDevID_3,PrimProtDevNum_3,TransformerID_3,TransformerNum_3,SubXs_3,SecProtDevID_3,SecProtDevNum_3,BranchStatus_3,BranchID_4,PrimProtDevID_4,PrimProtDevNum_4,TransformerID_4,TransformerNum_4,SubXs_4,SecProtDevID_4,SecProtDevNum_4,BranchStatus_4,BranchID_5,PrimProtDevID_5,PrimProtDevNum_5,TransformerID_5,TransformerNum_5,SubXs_5,SecProtDevID_5,SecProtDevNum_5,BranchStatus_5,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[TRANSFORMER]')) 
             ySUBSTATION=tline;
             fprintf(fid2,'%s\r\n',ySUBSTATION);
             tline = fgetl(fid);
       end
   end  
%    if ~isempty(strfind(tline,'FORMAT_SUBSTATION=ID,MVA,MVA_1,MVA_2,MVA_3,MVA_4,KVLL,KVLLdesired,R1,X1,R0,X0,Conn,PhaseAngle,PrimaryEquivalentType,SubEqVal1,SubEqVal2,SubEqVal3,SubEqVal4,SubPrimaryLLVoltage,SecondaryFaultReactance,TxfoConnection,HarmonicEnveloppe,ImpedanceUnit,BranchID_1,PrimProtDevID_1,PrimProtDevNum_1,TransformerID_1,TransformerNum_1,SubXs_1,SecProtDevID_1,SecProtDevNum_1,BranchStatus_1,BranchID_2,PrimProtDevID_2,PrimProtDevNum_2,TransformerID_2,TransformerNum_2,SubXs_2,SecProtDevID_2,SecProtDevNum_2,BranchStatus_2,BranchID_3,PrimProtDevID_3,PrimProtDevNum_3,TransformerID_3,TransformerNum_3,SubXs_3,SecProtDevID_3,SecProtDevNum_3,BranchStatus_3,BranchID_4,PrimProtDevID_4,PrimProtDevNum_4,TransformerID_4,TransformerNum_4,SubXs_4,SecProtDevID_4,SecProtDevNum_4,BranchStatus_4,BranchID_5,PrimProtDevID_5,PrimProtDevNum_5,TransformerID_5,TransformerNum_5,SubXs_5,SecProtDevID_5,SecProtDevNum_5,BranchStatus_5,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[TRANSFORMER]')) 
%              ySUBSTATION=tline;
%              fprintf(fid2,'%s\r\n',ySUBSTATION);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% TRANSFORMER data
fid2 = fopen([FileName '_TRANSFORMER.txt'], 'w');
if ~isempty(strfind(tline,'[TRANSFORMER]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_TRANSFORMER=ID,Type,KVA,KVLLprim,KVLLsec,Z1,Z0,XR,XR0,Conn,Rg_prim,Xg_prim,Rg_sec,Xg_sec,IsLTC,Taps,LowerBandwidth,UpperBandwidth,MinReg_Range,MaxReg_Range,Reversible,SelfCooledKVA,SelfCooledKVA_2,NoLoadLosses,SymbolID,PhaseShiftType,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[THREE WINDING TRANSFORMER]')) 
             yTRANSFORMER=tline;
             fprintf(fid2,'%s\r\n',yTRANSFORMER);
             tline = fgetl(fid);
       end
   end    
%    if ~isempty(strfind(tline,'FORMAT_TRANSFORMER=ID,Type,KVA,KVLLprim,KVLLsec,Z1,Z0,Z0PrimSec,Z0PrimMag,Z0SecMag,XR,XR0,XR0PrimSec,XR0PrimMag,XR0SecMag,Conn,Rg_prim,Xg_prim,Rg_sec,Xg_sec,IsLTC,Taps,LowerBandwidth,UpperBandwidth,MinReg_Range,MaxReg_Range,Reversible,SelfCooledKVA,SelfCooledKVA_2,SelfCooledKVA_3,SelfCooledKVA_4,NoLoadLosses,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,SymbolID,PhaseShiftType,Comments,DryType'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[THREE WINDING TRANSFORMER]')) 
%              yTRANSFORMER=tline;
%              fprintf(fid2,'%s\r\n',yTRANSFORMER);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% THREE WINDING TRANSFORMER data
fid2 = fopen([FileName '_THREEWINDINGTRANSFORMER.txt'], 'w');
if ~isempty(strfind(tline,'[THREE WINDING TRANSFORMER]'))
    tline = fgetl(fid);
    if ~isempty(strfind(tline,'FORMAT_THREEWINDINGTRANSFORMER=ID,PrimaryRatedCapacity,PrimaryCapacityLimit1,PrimaryCapacityLimit2,PrimaryVoltage,PrimaryConnection,PrimaryRg,PrimaryXg,PrimaryToSecondaryZ1,PrimaryToSecondaryZ0,PrimaryToSecondaryXR1,PrimaryToSecondaryXR0,PrimaryToSecondaryPhaseShiftType,PrimaryToTertiaryZ1,PrimaryToTertiaryZ0,PrimaryToTertiaryXR1,PrimaryToTertiaryXR0,PrimaryToTertiaryPhaseShiftType,SecondaryToTertiaryZ1,SecondaryToTertiaryZ0,SecondaryToTertiaryXR1,SecondaryToTertiaryXR0,SecondaryRatedCapacity,SecondaryCapacityLimit1,SecondaryCapacityLimit2,SecondaryVoltage,SecondaryConnection,SecondaryRg,SecondaryXg,TertiaryRatedCapacity,TertiaryCapacityLimit1,TertiaryCapacityLimit2,TertiaryVoltage,TertiaryConnection,TertiaryRg,TertiaryXg,LTC1_NumberOfTaps,LTC1_UpperBandwidth,LTC1_LowerBandwidth,LTC1_MaximumRegulationRange,LTC1_MinimumRegulationRange,LTC2_NumberOfTaps,LTC2_UpperBandwidth,LTC2_LowerBandwidth,LTC2_MaximumRegulationRange,LTC2_MinimumRegulationRange,SymbolID,Comments'))
        tline = fgetl(fid);
        while isempty(strfind(tline,'[SERIE CAPACITOR]'))
            yTHREEWINDINGTRANSFORMER = tline;
            fprintf(fid2,'%s\r\n',yTHREEWINDINGTRANSFORMER);
            tline = fgetl(fid);
        end
	end
%         if ~isempty(strfind(tline,'FORMAT_THREEWINDINGTRANSFORMER=ID,PrimaryRatedCapacity,PrimaryCapacityLimit1,PrimaryCapacityLimit2,PrimaryCapacityLimit3,PrimaryCapacityLimit4,PrimaryVoltage,PrimaryConnection,PrimaryRg,PrimaryXg,PrimaryToSecondaryZ1,PrimaryToSecondaryZ0,PrimaryToSecondaryXR1,PrimaryToSecondaryXR0,PrimaryToSecondaryPhaseShiftType,PrimaryToTertiaryZ1,PrimaryToTertiaryZ0,PrimaryToTertiaryXR1,PrimaryToTertiaryXR0,PrimaryToTertiaryPhaseShiftType,SecondaryToTertiaryZ1,SecondaryToTertiaryZ0,SecondaryToTertiaryXR1,SecondaryToTertiaryXR0,SecondaryRatedCapacity,SecondaryCapacityLimit1,SecondaryCapacityLimit2,SecondaryCapacityLimit3,SecondaryCapacityLimit4,SecondaryVoltage,SecondaryConnection,SecondaryRg,SecondaryXg,TertiaryRatedCapacity,TertiaryCapacityLimit1,TertiaryCapacityLimit2,TertiaryCapacityLimit3,TertiaryCapacityLimit4,TertiaryVoltage,TertiaryConnection,TertiaryRg,TertiaryXg,LTC1_NumberOfTaps,LTC1_UpperBandwidth,LTC1_LowerBandwidth,LTC1_MaximumRegulationRange,LTC1_MinimumRegulationRange,LTC2_NumberOfTaps,LTC2_UpperBandwidth,LTC2_LowerBandwidth,LTC2_MaximumRegulationRange,LTC2_MinimumRegulationRange,NoLoadLosses,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,SymbolID,Comments'))
%             tline = fgetl(fid);
%            while isempty(strfind(tline,'[AUTOTRANSFORMER]'))
%                yTHREEWINDINGTRANSFORMER = tline;
%                fprintf(fid2,'%s\r\n',yTHREEWINDINGTRANSFORMER);
%                tline = fgetl(fid);
%            end
%         end
   fclose(fid2);
end

% % AUTOTRANSFORMER data
% fid2 = fopen([FileName '_AUTOTRANSFORMER.txt'], 'w');
% if ~isempty(strfind(tline,'[AUTOTRANSFORMER]'))
%    tline = fgetl(fid);
%     if ~isempty(strfind(tline,'FORMAT_AUTOTRANSFORMER=ID,Type,KVA,KVLLprim,KVLLsec,Z1,Z0,Z0PrimSec,Z0PrimMag,Z0SecMag,XR,XR0,XR0PrimSec,XR0PrimMag,XR0SecMag,ConnectionConfiguration,Rg,Xg,IsLTC,Taps,LowerBandwidth,UpperBandwidth,MinReg_Range,MaxReg_Range,Reversible,SelfCooledKVA_1,SelfCooledKVA_2,SelfCooledKVA_3,SelfCooledKVA_4,NoLoadLosses,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,SymbolID,Comments'))
%        tline = fgetl(fid);
%            while isempty(strfind(tline,'[THREE WINDING AUTO TRANSFORMER]')) 
%                  yAUTOTRANSFORMER=tline;
%                  fprintf(fid2,'%s\r\n',yAUTOTRANSFORMER);
%                  tline = fgetl(fid);
%            end   
%    end    
%    fclose(fid2);
% end
% % THREE WINDING AUTO TRANSFORMER data
% fid2 = fopen([FileName '_THREEWINDINGAUTOTRANSFORMER.txt'], 'w');
% if ~isempty(strfind(tline,'[THREE WINDING AUTO TRANSFORMER]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_THREEWINDINGAUTOTRANSFORMER=ID,PrimaryRatedCapacity,PrimaryCapacityLimit1,PrimaryCapacityLimit2,PrimaryCapacityLimit3,PrimaryCapacityLimit4,PrimaryVoltage,PrimarySecondaryConnection,PrimaryRg,PrimaryXg,PrimaryToSecondaryZ1,PrimaryToSecondaryZ0,PrimaryToSecondaryXR1,PrimaryToSecondaryXR0,PrimaryToTertiaryZ1,PrimaryToTertiaryZ0,PrimaryToTertiaryXR1,PrimaryToTertiaryXR0,SecondaryToTertiaryZ1,SecondaryToTertiaryZ0,SecondaryToTertiaryXR1,SecondaryToTertiaryXR0,SecondaryRatedCapacity,SecondaryCapacityLimit1,SecondaryCapacityLimit2,SecondaryCapacityLimit3,SecondaryCapacityLimit4,SecondaryVoltage,SecondaryRg,SecondaryXg,TertiaryRatedCapacity,TertiaryCapacityLimit1,TertiaryCapacityLimit2,TertiaryCapacityLimit3,TertiaryCapacityLimit4,TertiaryVoltage,TertiaryConnection,TertiaryRg,TertiaryXg,LTC1_NumberOfTaps,LTC1_UpperBandwidth,LTC1_LowerBandwidth,LTC1_MaximumRegulationRange,LTC1_MinimumRegulationRange,LTC2_NumberOfTaps,LTC2_UpperBandwidth,LTC2_LowerBandwidth,LTC2_MaximumRegulationRange,LTC2_MinimumRegulationRange,NoLoadLosses,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,SymbolID,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SERIE CAPACITOR]')) 
%              yTHREEWINDINGAUTOTRANSFORMER=tline;
%              fprintf(fid2,'%s\r\n',yTHREEWINDINGAUTOTRANSFORMER);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end

% SERIE CAPACITOR data
fid2 = fopen([FileName '_SERIECAPACITOR.txt'], 'w');
if ~isempty(strfind(tline,'[SERIE CAPACITOR]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SERIECAPACITOR=ID,Amps,Amps_2,Reactance,SymbolID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SERIE REACTOR]')) 
             ySERIECAPACITOR=tline;
             fprintf(fid2,'%s\r\n',ySERIECAPACITOR);
             tline = fgetl(fid);
       end
   end  
%    if ~isempty(strfind(tline,'FORMAT_SERIECAPACITOR=ID,Amps,Amps_1,Amps_2,Amps_3,Amps_4,Reactance,SymbolID,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SERIE REACTOR]')) 
%              ySERIECAPACITOR=tline;
%              fprintf(fid2,'%s\r\n',ySERIECAPACITOR);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% SERIE REACTOR data
fid2 = fopen([FileName '_SERIEREACTOR.txt'], 'w');
if ~isempty(strfind(tline,'[SERIE REACTOR]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SERIEREACTOR=ID,Amps,Amps_2,Reactance,SymbolID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SHUNT CAPACITOR]')) 
             ySERIEREACTOR=tline;
             fprintf(fid2,'%s\r\n',ySERIEREACTOR);
             tline = fgetl(fid);
       end
   end    
%    if ~isempty(strfind(tline,'FORMAT_SERIEREACTOR=ID,Amps,Amps_1,Amps_2,Amps_3,Amps_4,Reactance,SymbolID,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SHUNT CAPACITOR]')) 
%              ySERIEREACTOR=tline;
%              fprintf(fid2,'%s\r\n',ySERIEREACTOR);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% SHUNT CAPACITOR data
fid2 = fopen([FileName '_SHUNTCAPACITOR.txt'], 'w');
if ~isempty(strfind(tline,'[SHUNT CAPACITOR]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SHUNTCAPACITOR=ID,KVAR,KV,CostForFixedBank,CostForSwitchedBank,SymbolID,Type,InterruptingRating,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SHUNT REACTOR]')) 
             ySHUNTCAPACITOR=tline;
             fprintf(fid2,'%s\r\n',ySHUNTCAPACITOR);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% SHUNT REACTOR data
fid2 = fopen([FileName '_SHUNTREACTOR.txt'], 'w');
if ~isempty(strfind(tline,'[SHUNT REACTOR]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SHUNTREACTOR=ID,KVAR,KV,CostForFixedBank,CostForSwitchedBank,SymbolID,Type,InterruptingRating,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SHUNT CAPACITOR]')) 
             ySHUNTREACTOR=tline;
             fprintf(fid2,'%s\r\n',ySHUNTREACTOR);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% SHUNT CAPACITOR data
fid2 = fopen([FileName '_SHUNTCAPACITOR.txt'], 'w');
if ~isempty(strfind(tline,'[SHUNT CAPACITOR]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SHUNTCAPACITOR=ID,KVAR,KVLN,CostForFixedBank,CostForSwitchedBank,SymbolID,Type,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[ARC FURNACE]')) 
             ySHUNTCAPACITOR=tline;
             fprintf(fid2,'%s\r\n',ySHUNTCAPACITOR);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% ARC FURNACE data
fid2 = fopen([FileName '_ARCFURNACE.txt'], 'w');
if ~isempty(strfind(tline,'[ARC FURNACE]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_ARCFURNACE=ID,Frequency,MagnitudePercent,SymbolOnID,SymbolOffID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[CTYPE FILTER]')) 
             yARCFURNACE=tline;
             fprintf(fid2,'%s\r\n',yARCFURNACE);
             tline = fgetl(fid);
       end
   end    
   
   fclose(fid2);
end

% CTYPE FILTER data
fid2 = fopen([FileName '_CTYPEFILTER.txt'], 'w');
if ~isempty(strfind(tline,'[CTYPE FILTER]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_CTYPEFILTER=ID,R1,L1,C1,R2,L2,C2,R3,L3,C3,SymbolOnID,SymbolOffID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[DOUBLE TUNED FILTER]')) 
             yCTYPEFILTER=tline;
             fprintf(fid2,'%s\r\n',yCTYPEFILTER);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% DOUBLE TUNED FILTER data
fid2 = fopen([FileName '_DOUBLETUNEDFILTER.txt'], 'w');
if ~isempty(strfind(tline,'[DOUBLE TUNED FILTER]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_DOUBLETUNEDFILTER=ID,R1,L1,C1,R2,L2,C2,R3,TunedFrequency1,TunedFrequency2,FirstTunedQ,FirstTunedQualityFactor,FirstTunedV,SecondTunedQ,SecondTunedQualityFactor,SecondTunedV,SymbolOnID,SymbolOffID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[HIGH PASS FILTER]')) 
             yDOUBLETUNEDFILTER=tline;
             fprintf(fid2,'%s\r\n',yDOUBLETUNEDFILTER);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% HIGH PASS FILTER data
fid2 = fopen([FileName '_HIGHPASSFILTER.txt'], 'w');
if ~isempty(strfind(tline,'[HIGH PASS FILTER]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_HIGHPASSFILTER=ID,R,L,C,KVAR,KV,SymbolOnID,SymbolOffID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[IDEAL CONVERTER]')) 
             yHIGHPASSFILTER=tline;
             fprintf(fid2,'%s\r\n',yHIGHPASSFILTER);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% IDEAL CONVERTER data
fid2 = fopen([FileName '_IDEALCONVERTER.txt'], 'w');
if ~isempty(strfind(tline,'[IDEAL CONVERTER]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_IDEALCONVERTER=ID,PulseNumber,KVA,KV,SymbolOnID,SymbolOffID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[NON IDEAL CONVERTER]')) 
             yIDEALCONVERTER=tline;
             fprintf(fid2,'%s\r\n',yIDEALCONVERTER);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% NON IDEAL CONVERTER data
fid2 = fopen([FileName '_NONIDEALCONVERTER.txt'], 'w');
if ~isempty(strfind(tline,'[NON IDEAL CONVERTER]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_NONIDEALCONVERTER=ID,PulseNumber,KVA,KV,SymbolOnID,SymbolOffID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SHUNT MULTI FREQUENCY SOURCE]')) 
             yNONIDEALCONVERTER=tline;
             fprintf(fid2,'%s\r\n',yNONIDEALCONVERTER);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% SHUNT MULTI FREQUENCY SOURCE data
fid2 = fopen([FileName '_SHUNTMULTIFREQUENCYSOURCE.txt'], 'w');
if ~isempty(strfind(tline,'[SHUNT MULTI FREQUENCY SOURCE]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SHUNTMULTIFREQUENCYSOURCE=ID,Frequency,Amp,Angle,InAmp,FundamentalCurrent,SymbolOnID,SymbolOffID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SHUNT SINGLE FREQUENCY SOURCE]')) 
             ySHUNTMULTIFREQUENCYSOURCE=tline;
             fprintf(fid2,'%s\r\n',ySHUNTMULTIFREQUENCYSOURCE);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% SHUNT MULTI FREQUENCY SOURCE data
fid2 = fopen([FileName '_SHUNTSINGLEFREQUENCYSOURCE.txt'], 'w');
if ~isempty(strfind(tline,'[SHUNT SINGLE FREQUENCY SOURCE]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SHUNTSINGLEFREQUENCYSOURCE=ID,Frequency,Amp,Angle,SymbolOnID,SymbolOffID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SINGLE TUNED FILTER]')) 
             ySHUNTSINGLEFREQUENCYSOURCE=tline;
             fprintf(fid2,'%s\r\n',ySHUNTSINGLEFREQUENCYSOURCE);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% % FREQUENCY SOURCE data
% fid2 = fopen([FileName '_FREQUENCYSOURCE.txt'], 'w');
% if ~isempty(strfind(tline,'[FREQUENCY SOURCE]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_FREQUENCYSOURCE=ID,Frequency,Magnitude,Angle,InPercent,SymbolOnID,SymbolOffID,VoltageSource,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SINGLE TUNED FILTER]')) 
%              yFREQUENCYSOURCE=tline;
%              fprintf(fid2,'%s\r\n',yFREQUENCYSOURCE);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end

% SINGLE TUNED FILTER data
fid2 = fopen([FileName '_SINGLETUNEDFILTER.txt'], 'w');
if ~isempty(strfind(tline,'[SINGLE TUNED FILTER]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SINGLETUNEDFILTER=ID,R,L,C,TunedFrequency,KVAR,KV,ConnectionConfiguration,QualityFactor,SymbolOnID,SymbolOffID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[GROUNDING TRANSFORMER]')) 
             ySINGLETUNEDFILTER=tline;
             fprintf(fid2,'%s\r\n',ySINGLETUNEDFILTER);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% GROUNDING TRANSFORMER data
fid2 = fopen([FileName '_GROUNDINGTRANSFORMER.txt'], 'w');
if ~isempty(strfind(tline,'[GROUNDING TRANSFORMER]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_GROUNDINGTRANSFORMER=ID,Z1,Z0,X1_R1,X0_R0,RatedCapacity,RatedVoltage,ConnectionConfiguration,SymbolID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[WECS]')) 
             yGROUNDINGTRANSFORMER=tline;
             fprintf(fid2,'%s\r\n',yGROUNDINGTRANSFORMER);
             tline = fgetl(fid);
       end
   end  
%    if ~isempty(strfind(tline,'FORMAT_GROUNDINGTRANSFORMER=ID,Z1,Z0,X1_R1,X0_R0,RatedCapacity,RatedVoltage,ConnectionConfiguration,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,SymbolID,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[WECS]')) 
%              yGROUNDINGTRANSFORMER=tline;
%              fprintf(fid2,'%s\r\n',yGROUNDINGTRANSFORMER);
%              tline = fgetl(fid);
%        end
%    end    
   fclose(fid2);
end

% WECS data
fid2 = fopen([FileName '_WECS.txt'], 'w');
if ~isempty(strfind(tline,'[WECS]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_WECS=ID,RatedPower_KW,MaxPower_KW,RatedWindspeed_MS,CutinWindspeed_MS,CutoutWindspeed_MS,NbOfRotorBlades,ROTORRADIUS,RatedSpeed_RPM,MinSpeed_RPM,MaxSpeed_RPM,WindTurbineInertia_KGM2,SpringConstant_NMSRAD,DampingConstant_NMSRAD,GearboxRatio,GeneratorType,GeneratorRatedCapacity_KVA,GeneratorRatedVoltage_KV,GeneratorRatedPower_KW,GeneratorPowerFactor_Percent,GeneratorEfficiency_Percent,GeneratorRatedSpeed_RPM,SymbolID,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[INDUCTION MACHINE EQ CIRCUIT]')) 
             yWECS=tline;
             fprintf(fid2,'%s\r\n',yWECS);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% INDUCTION MACHINE EQ CIRCUIT data
fid2 = fopen([FileName '_INDUCTIONMACHINEEQCIRCUIT.txt'], 'w');
if ~isempty(strfind(tline,'[INDUCTION MACHINE EQ CIRCUIT]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_INDUCTIONMACHINEEQCIRCUIT=ID,Type,RotorType,EstimationMethod,StatorRS_Ohms,StatorXS_Ohms,MagnetisingRM_Ohms,MagnetisingXM_Ohms,OuterCageRotorRR1_OHMS,OuterCageRotorXR1_OHMS,InnerCageRotorRR2_OHMS,InnerCageRotorXR2_OHMS,CageFactorR,CageFactorX,InertiaUnit,Inertia,ImpedanceUnit'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SYNCHRONOUS MACHINE EQ CIRCUIT]')) 
             yINDUCTIONMACHINEEQCIRCUIT=tline;
             fprintf(fid2,'%s\r\n',yINDUCTIONMACHINEEQCIRCUIT);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% SYNCHRONOUS MACHINE EQ CIRCUIT data
fid2 = fopen([FileName '_SYNCHRONOUSMACHINEEQCIRCUIT.txt'], 'w');
if ~isempty(strfind(tline,'[SYNCHRONOUS MACHINE EQ CIRCUIT]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SYNCHRONOUSMACHINEEQCIRCUIT=ID,Type,SynchronousXD,SynchronousXQ,SynchronousXP,TransientXD,TransientXQ,TransientTDO,TransientTQO,SubtransientXD,SubtransientXQ,SubtransientTDO,SubtransientTQO,InertiaUnits,Inertia'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SYNCHRONOUS MACHINE EXT STAB]')) 
         ySYNCHRONOUSMACHINEEQCIRCUIT=tline;
         fprintf(fid2,'%s\r\n',ySYNCHRONOUSMACHINEEQCIRCUIT);
         tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end
 
% SYNCHRONOUS MACHINE EXT STAB data
fid2 = fopen([FileName '_SYNCHRONOUSMACHINEEXTSTAB.txt'], 'w');
if ~isempty(strfind(tline,'[SYNCHRONOUS MACHINE EXT STAB]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SYNCHRONOUSMACHINEEXTST=ID,Type,StabilityModel,DampingConstant,SaturationSGU,SaturationSGL,SaturationEU,SaturationEL'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[SYNCHRONOUS MACHINE EXT HARMO]')) 
             ySYNCHRONOUSMACHINEEXTSTAB=tline;
             fprintf(fid2,'%s\r\n',ySYNCHRONOUSMACHINEEXTSTAB);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% SYNCHRONOUS MACHINE EXT HARMO data
fid2 = fopen([FileName '_SYNCHRONOUSMACHINEEXTHARMO.txt'], 'w');
if ~isempty(strfind(tline,'[SYNCHRONOUS MACHINE EXT HARMO]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_SYNCHRONOUSMACHINEEXTHA=ID,Type,ComputeMode,R2ohms,L2mh,SkinEffectA,SkinEffectB'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[WIND MODEL]')) 
             ySYNCHRONOUSMACHINEEXTHARMO=tline;
             fprintf(fid2,'%s\r\n',ySYNCHRONOUSMACHINEEXTHARMO);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% WIND MODEL data
fid2 = fopen([FileName '_WINDMODEL.txt'], 'w');
if ~isempty(strfind(tline,'[WIND MODEL]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_WINDMODEL=ID,FromFile,FileName,Comments'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[WIND MODEL POINTS]')) 
             yWINDMODEL=tline;
             fprintf(fid2,'%s\r\n',yWINDMODEL);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% WIND MODEL POINTS data
fid2 = fopen([FileName '_WINDMODELPOINTS.txt'], 'w');
if ~isempty(strfind(tline,'[WIND MODEL POINTS]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_WINDMODELPOINTS=ID,PointIndex,Time_S,WindSpeed_MS'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'   ')) 
             yWINDMODEL=tline;
             fprintf(fid2,'%s\r\n',yWINDMODEL);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% % INSOLATION MODEL data
% fid2 = fopen([FileName '_INSOLATIONMODEL.txt'], 'w');
% if ~isempty(strfind(tline,'[INSOLATION MODEL]'))
% 
%    tline = fgetl(fid);
%     if ~isempty(strfind(tline,'FORMAT_INSOLATIONMODEL=ID,FromFile,FileName,Comments'))
% 
%         tline = fgetl(fid);
%        while isempty(strfind(tline,'[PHOTOVOLTAIC]'))
% 
%            yINSOLATIONMODEL = tline;
%            fprintf(fid2,'%s\r\n',yINSOLATIONMODEL);
%            tline = fgetl(fid);
%        end
%     end
%    fclose(fid2);
% end
% 
% % PHOTOVOLTAIC data
% fid2 = fopen([FileName '_PHOTOVOLTAIC.txt'], 'w');
% if ~isempty(strfind(tline,'[PHOTOVOLTAIC]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_PHOTOVOLTAIC=ID,MPPCurrent,MPPVoltage,SCCurrent,OCVoltage,SCCurrentTempCoeff,OCVoltageTempCoeff,OperatingTemperature,RefAmbientTemperature,STCTemperature,STCInsolation,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SOFC]')) 
%              yPHOTOVOLTAIC=tline;
%              fprintf(fid2,'%s\r\n',yPHOTOVOLTAIC);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% % SOFC data
% fid2 = fopen([FileName '_SOFC.txt'], 'w');
% if ~isempty(strfind(tline,'[SOFC]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_SOFC=ID,RatedPower,N0,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[MICRO TURBINE]')) 
%              ySOFC=tline;
%              fprintf(fid2,'%s\r\n',ySOFC);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% % MICRO TURBINE data
% fid2 = fopen([FileName '_MICROTURBINE.txt'], 'w');
% if ~isempty(strfind(tline,'[MICRO TURBINE]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_MICROTURBINE=ID,GovernorKp,GovernorKi,TurbineTimeConstant,Inertia,InertiaUnitType,RatedCapacity,RatedVoltage,RatedPower,RatedSpeed,SynchronousReactance,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[NETWORKPROTECTOR]')) 
%              yMICROTURBINE=tline;
%              fprintf(fid2,'%s\r\n',yMICROTURBINE);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% % NETWORKPROTECTOR data
% fid2 = fopen([FileName '_NETWORKPROTECTOR.txt'], 'w');
% if ~isempty(strfind(tline,'[NETWORKPROTECTOR]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_NETWORKPROTECTOR=ID,Amps,Amps_1,Amps_2,Amps_3,Amps_4,KVLL,Reversible,InterruptingRating,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,StuckProbability,SwitchTime,SymbolOpenID,SymbolCloseID,SinglePhaseLocking,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[SVC]')) 
%              yNETWORKPROTECTOR=tline;
%              fprintf(fid2,'%s\r\n',yNETWORKPROTECTOR);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% 
% % SVC data
% fid2 = fopen([FileName '_SVC.txt'], 'w');
% if ~isempty(strfind(tline,'[SVC]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_SVC=ID,PulseNumber,RatedVoltage,MaxReactivePower,MinReactivePower,SymbolOnID,SymbolOffID,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[GENERATOR COST CURVE MODEL]')) 
%              ySVC=tline;
%              fprintf(fid2,'%s\r\n',ySVC);
%              tline = fgetl(fid);
%              if ~ischar(tline),   break
%              end
%        end
%    end    
%    fclose(fid2);
% end
% 
% % GENERATOR COST CURVE MODEL data
% fid2 = fopen([FileName '_GENERATORCOSTCURVEMODEL.txt'], 'w');
% if ~isempty(strfind(tline,'[GENERATOR COST CURVE MODEL]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_GENERATORCOSTCURVEMODEL=ID,SymbolOnID,SymbolOffID,Favorite,GeneratorCostCurveType,IntegrationConstant,LinearCoefficientPerKW,QuadraticCoefficientPerKW2,ExponentialCoefficient,ExponentialScaleFactorPerKW,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[BUSWAY]')) 
%              yGENERATORCOSTCURVEMODEL=tline;
%              fprintf(fid2,'%s\r\n',yGENERATORCOSTCURVEMODEL);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% % BUSWAY data
% fid2 = fopen([FileName '_BUSWAY.txt'], 'w');
% if ~isempty(strfind(tline,'[BUSWAY]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_BUSWAY=ID,SymbolOnID,SymbolOffID,Favorite,RatedVoltage,NominalRating,FirstRating,SecondRating,ThirdRating,FourthRating,WithstandRating,Material,BuswayType,PositiveSequenceResistance,PositiveSequenceReactance,ZeroSequenceResistance,ZeroSequenceReactance,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[PHASE SHIFTER TRANSFORMER]')) 
%              yBUSWAY=tline;
%              fprintf(fid2,'%s\r\n',yBUSWAY);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% % GENERATOR COST CURVE MODEL data
% fid2 = fopen([FileName '_PHASESHIFTERTRANSFORMER.txt'], 'w');
% if ~isempty(strfind(tline,'[PHASE SHIFTER TRANSFORMER]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_PHASESHIFTERTRANSFORMER=ID,SymbolOnID,SymbolOffID,Favorite,Reversible,RatedVoltage,NominalRating,FirstLoadingLimit,SecondLoadingLimit,ThirdLoadingLimit,FourthLoadingLimit,PositiveSequenceImpedance,ZeroSequenceImpedancePrimary,ZeroSequenceImpedanceSecondary,ZeroSequenceImpedanceTertiary,XRRatio,XR0RatioPrimary,XR0RatioSecondary,XR0RatioTertiary,MinimumPhaseShift,MaximumPhaseShift,NumberOfTaps,FailRate,TmpFailRate,MajorRepairTime,MinorRepairTime,MajorFailureProportion,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[GENERATION CURVE MODEL]')) 
%              yPHASESHIFTERTRANSFORMER=tline;
%              fprintf(fid2,'%s\r\n',yPHASESHIFTERTRANSFORMER);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% % GENERATION CURVE MODEL data
% fid2 = fopen([FileName '_GENERATIONCURVEMODEL.txt'], 'w');
% if ~isempty(strfind(tline,'[GENERATION CURVE MODEL]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_GENERATIONCURVEMODEL=ID,SymbolOnID,SymbolOffID,Favorite,CurveModel,FromFile,FileName,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[LOAD CURVE MODEL]')) 
%              yGENERATIONCURVEMODEL=tline;
%              fprintf(fid2,'%s\r\n',yGENERATIONCURVEMODEL);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% % LOAD CURVE MODEL data
% fid2 = fopen([FileName '_LOADCURVEMODEL.txt'], 'w');
% if ~isempty(strfind(tline,'[LOAD CURVE MODEL]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_LOADCURVEMODEL=ID,SymbolOnID,SymbolOffID,Favorite,CurveModel,FromFile,FileName,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[MOTOR CURVE MODEL]')) 
%              yLOADCURVEMODEL=tline;
%              fprintf(fid2,'%s\r\n',yLOADCURVEMODEL);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% % GENERATOR COST CURVE MODEL data
% fid2 = fopen([FileName '_MOTORCURVEMODEL.txt'], 'w');
% if ~isempty(strfind(tline,'[MOTOR CURVE MODEL]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_MOTORCURVEMODEL=ID,SymbolOnID,SymbolOffID,Favorite,CurveModel,FromFile,FileName,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'[VARIABLE FREQUENCY DRIVE]')) 
%              yMOTORCURVEMODEL=tline;
%              fprintf(fid2,'%s\r\n',yMOTORCURVEMODEL);
%              tline = fgetl(fid);
%        end
%    end    
%    fclose(fid2);
% end
% 
% % VARIABLE FREQUENCY DRIVE data
% fid2 = fopen([FileName '_VARIABLEFREQUENCYDRIVE.txt'], 'w');
% if ~isempty(strfind(tline,'[VARIABLE FREQUENCY DRIVE]'))
%    tline = fgetl(fid);
%    if ~isempty(strfind(tline,'FORMAT_VARIABLEFREQUENCYDRIVE=ID,Amps,Amps_1,Amps_2,Amps_3,Amps_4,KVLL,RatedPower,Efficiency,SymbolOpenID,SymbolCloseID,Favorite,Comments'))
%        tline = fgetl(fid);
%        while isempty(strfind(tline,'    ')) 
%              yVARIABLEFREQUENCYDRIVE=tline;
%              fprintf(fid2,'%s\r\n',yVARIABLEFREQUENCYDRIVE);
%              tline = fgetl(fid);
%              if ~ischar(tline),   break
%              end
%        end
%    end    
%    fclose(fid2);
% end



ExecutionFlag = 1;

fclose(fid);



    