function[struc] = GetCYMData_transformersetting(FilePath,flgTest)

if nargin<2
    flgTest=0;
end

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
% 1SectionID,2Location,3EqID,4DeviceNumber,5CoordX,6CoordY,7Conn,8PrimTap,
% 9SecondaryTap,10RgPrim,11XgPrim,12RgSec,13XgSec,14ODPrimPh,15PrimaryBaseVoltage,
% 16SecondaryBaseVoltage,17FromNodeID,18SettingOption,19SetPoint,20ControlType,
% 21LowerBandWidth,22UpperBandWidth,23TapLocation,24InitialTapPosition,25InitialTapPositionMode,
% 26Tap,27MaxBuck,28MaxBoost,29CT,30PT,31Rset,32Xset,33FirstHouseHigh,34FirstHouseLow,35PhaseON,
% 36AtSectionID, 37MasterID,38FaultIndicator,39DiversityFactorA,40DiversityFactorB,41DiversityFactorC,
% 42PhaseShiftType,43GammaPhaseShift,44DemandType,45Value1A,46Value2A,47Value1B,48Value2B,
% 49Value1C,50Value2C,51Value1ABC,52Value2ABC,53ConnectionStatus,54Reversible
strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
if flgTest
    FilePath=[FileName '_TRANSFORMERSETTING_small.txt'];
else
    FilePath=[FileName '_TRANSFORMERSETTING.txt'];
end

if ~exist(FilePath)
    disp('"transformersetting" structure not populated');
    disp('File not found.');
    FilePath
    struc=[];
    return
else
    fid = fopen(FilePath, 'r');
    y=textscan(fid,strFormat, 'CollectOutput', 0,'delimiter',strDelimiter,'MultipleDelimsAsOne',0);
    fclose(fid);
end
struc.SectionID=fnSanitize(y{1});
struc.Location=y{2};
struc.EqID=y{3};
struc.DeviceNumber=y{4};
struc.Conn=y{7};
struc.tapVoltageprim=y{8};
struc.tapVoltagesec=y{9};
struc.Rgprim=y{10};
struc.Xgprim=y{11};
struc.Rgsec=y{12};
struc.Xgsec=y{13};






    