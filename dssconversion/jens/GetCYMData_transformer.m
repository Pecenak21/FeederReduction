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
% Format:
% 1ID,2Type,3KVA,4KVLLprim,5KVLLsec,6Z1,7Z0,8Z0PrimSec,9Z0PrimMag,10Z0SecMag,
% 11XR,12XR0,13XR0PrimSec,14XR0PrimMag,15XR0SecMag,16Conn,17Rg_prim,18Xg_prim,
% 19Rg_sec,20Xg_sec,21IsLTC,22Taps,23LowerBandwidth,24UpperBandwidth,25MinReg_Range,
% 26MaxReg_Range,27Reversible,28SelfCooledKVA,29SelfCooledKVA_2,30SelfCooledKVA_3,
% 31SelfCooledKVA_4,32NoLoadLosses,33FailRate,34TmpFailRate,35MajorRepairTime,
% 36MinorRepairTime,37MajorFailureProportion,38SymbolID,39PhaseShiftType,40Comments,41DryType

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_TRANSFORMER.txt'];

if ~exist(FilePath)
    disp('"tranformer" structure not populated');
    disp('File not found.');
    FilePath
    struc=[];
    return
else
    fid = fopen(FilePath, 'r');
    y=textscan(fid,strFormat, 'CollectOutput', 0,'delimiter',strDelimiter,'MultipleDelimsAsOne',0);
    fclose(fid);
end
struc.EqID=y{1};
struc.Type=y{2};
struc.KVA=y{3};
struc.KVLLprim=y{4};
struc.KVLLsec=y{5};
struc.Z1=y{6};
struc.Z0=y{7};
struc.XR=y{11};
struc.XR0=y{12};
struc.XR0PrimSec=y{13};
struc.XR0PrimMag=y{14};
struc.XR0SecMag=y{14};
struc.Conn=y{16};






    