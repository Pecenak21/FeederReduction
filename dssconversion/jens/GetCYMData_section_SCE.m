function[struc_Section, struc_Feeder] = GetCYMData_section_SCE(FilePath,flgTest)
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
% FORMAT_SECTION=1SectionID,2FromNodeID,3ToNodeID,4Phase,5ZoneID,6SubNetworkId,7EnvironmentID
% FORMAT_FEEDER=1NetworkID,2HeadNodeID,3CoordSet,4Year,5Description,6Color,
% 7LoadFactor,8Group1,9Group2,10TagText,11TagProperties,12TagDeltaX,
% 13TagDeltaY,14TagAngle,15TagAlignment,16TagBorder,17TagBackground,18TagTextColor,
% 19TagBorderColor,20TagBackgroundColor,21TagLocation,22TagFont,23TagTextSize,24Version,25EnvironmentID

% HydroOttawa format:
% FORMAT_FEEDER=1NetworkID,2HeadNodeID,3CoordSet,4Year,5Description,6Color,
% 7LoadFactor,8Group1,9Group2,10Group3,11TagText,12TagProperties,13TagDeltaX,
% 14TagDeltaY,15TagAngle,16TagAlignment,17TagBorder,18TagBackground,19TagTextColor,
% 20TagBorderColor,21TagBackgroundColor,22TagLocation,23TagFont,24TagTextSize,
% 25TagOffset,26Version,27EnvironmentID

strFormat_Section='%s %s %s %s %s %s %s';
strFormat_Feeder='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s';
strDelimiter=',';

FileName=RemoveFileExtension(FilePath);

if flgTest
    FilePath_Section=[FileName '_SECTION_small.txt'];
else
    FilePath_Section=[FileName '_SECTION.txt'];
end

FilePath_Feeder=[FileName '_FEEDER.txt'];

% Processing section data
if ~exist(FilePath_Section)
    disp('"Section" structure not populated.');
    disp('File not found.');
    FilePath
    struc=[];
    return
else
    fid = fopen(FilePath_Section, 'r');
    y=textscan(fid,strFormat_Section, 'CollectOutput', 0,'CommentStyle','FEEDER=','delimiter',strDelimiter,'MultipleDelimsAsOne',0);
    fclose(fid);
end
struc_Section.SectionID=fnSanitize(y{1});
struc_Section.FromNode=fnSanitize(y{2});
struc_Section.ToNode=fnSanitize(y{3});
struc_Section.Phasing=y{4};
struc_Section.NeutIsGrounded=0; %not sure if CYME provides this information (SynerGEE does), assuming that all Neutrals are not grounded

% does not read feeder data, might not be needed
% struc_Feeder.dummy=0;

% Processing feeder data
if ~exist(FilePath_Section)
    str=[FilePath_Section ' does not exist.'];
    return
else
    fid = fopen(FilePath_Section, 'r');
    y=textscan(fid,strFormat_Feeder, 'CollectOutput', 0,'delimiter',strDelimiter,'MultipleDelimsAsOne',0);
    fclose(fid);
end
strmatch('FEEDER=',y{1});

struc_Feeder.NetworkID=y{1};
struc_Feeder.HeadNodeID=y{2};
struc_Feeder.LoadFactor=y{7};
struc_Feeder.SubNetworkID=y{8};
struc_Feeder.kV=y{9};

index=strmatch('FEEDER=',struc_Feeder.NetworkID);
struc_Feeder=structreduce(struc_Feeder,index);






    