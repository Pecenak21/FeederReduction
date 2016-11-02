function[struc] = GetCYMData_transformer(FilePath)


% GetPSSEData_Bus(FilePath)
%
% PURPOSE : get content of bus parameters in file
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
%
% [SUBSTATION]
% FORMAT_SUBSTATION=1ID,2MVA,3MVA_1,4MVA_2,5MVA_3,6MVA_4,7KVLL,8KVLLdesired,9R1,10X1,
% 11R0,12X0,13Conn,14PhaseAngle,15PrimaryEquivalentType,16SubEqVal1,17SubEqVal2,
% 18SubEqVal3,19SubEqVal4,20SubPrimaryLLVoltage,21SecondaryFaultReactance,
% 22TxfoConnection,23HarmonicEnveloppe,24ImpedanceUnit,25BranchID_1,
% 26PrimProtDevID_1,27PrimProtDevNum_1,28TransformerID_1,29TransformerNum_1,
% 30SubXs_1,31SecProtDevID_1,32SecProtDevNum_1,33BranchStatus_1,34BranchID_2,
% 35PrimProtDevID_2,36PrimProtDevNum_2,37TransformerID_2,38TransformerNum_2,
% 39SubXs_2,40SecProtDevID_2,41SecProtDevNum_2,42BranchStatus_2,43BranchID_3,
% 44PrimProtDevID_3,45PrimProtDevNum_3,46TransformerID_3,47TransformerNum_3,
% 48SubXs_3,49SecProtDevID_3,50SecProtDevNum_3,51BranchStatus_3,52BranchID_4,
% 53PrimProtDevID_4,54PrimProtDevNum_4,55TransformerID_4,56TransformerNum_4,
% 57SubXs_4,58SecProtDevID_4,59SecProtDevNum_4,60BranchStatus_4,61BranchID_5,
% 62PrimProtDevID_5,63PrimProtDevNum_5,64TransformerID_5,65TransformerNum_5,
% 66SubXs_5,67SecProtDevID_5,68SecProtDevNum_5,69BranchStatus_5,70FailRate,
% 71TmpFailRate,72MajorRepairTime,73MinorRepairTime,74MajorFailureProportion,75Comments

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_SUBSTATION.txt'];

if ~exist(FilePath)
    disp('"substation" structure not populated');
    disp('File not found.');
    FilePath
    struc=[];
    return
else
    fid = fopen(FilePath, 'r');
    y=textscan(fid,strFormat, 'CollectOutput', 0,'delimiter',strDelimiter,'MultipleDelimsAsOne',0);
    fclose(fid);
end
struc.ID=fnSanitize(y{1});
struc.MVA=y{2};
struc.KVLL=y{7};
struc.R1=y{9};
struc.X1=y{10};
struc.R0=y{11};
struc.X0=y{12};
struc.Conn=y{13};
struc.PhaseAngle=y{14};
struc.TxfoConnection=y{22};
struc.ImpedanceUnit=y{24};







    