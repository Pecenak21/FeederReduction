function[struc] = GetCYMData_customerclass(FilePath)


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
% 1ID,2Description,3Color,4ConstantPower,5ConstantCurrent,6ConstantImpedance,
% 7UtilizationFactor,8PowerFactor,9MeteredLoads,10NonMeteredLoads,
% 11LoadFactor,12FrequencySensitivityP,13FrequencySensitivityQ,
% 14EnableHarmonic,15IsExponentialModel,16ConstantImpedanceZP,17ConstantImpedanceZQ,
% 18ConstantCurrentIP,19ConstantCurrentIQ,20ConstantPowerPP,
% 21ConstantPowerPQ,22ExponentialModelP,23ExponentialModelQ,
% 24LoadFlowVoltagePercentOfNominal,25AdjustmentSettings,26PowerCurveModel,27PowerCurveModelId 

strFormat='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);
FilePath=[FileName '_CUSTOMERCLASS.txt'];

if ~exist(FilePath)
    disp('"customerclass" structure not populated.');
    disp('File not found.');
    FilePath
    struc=[];
    return
else
    fid = fopen(FilePath, 'r');
    y=textscan(fid,strFormat, 'CollectOutput', 0,'delimiter',strDelimiter,'MultipleDelimsAsOne',0);
    fclose(fid);
end
struc.CustomerType=y{1};
struc.ConstantPower=y{4};
struc.ConstantCurrent=y{5};
struc.ConstantImpedance=y{6};
struc.UtilizationFactor=y{7};
struc.PowerFactor=y{8};
struc.LoadFactor=y{11};
struc.LoadFlowVoltagePercentOfNominal=y{24};






    