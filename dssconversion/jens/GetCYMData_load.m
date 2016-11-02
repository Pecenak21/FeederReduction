function[ExecutionFlag] = GetCYMData_load(FilePath)


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


yCUSTOMERCLASS=[];
yLOADS=[];
yCUSTOMERLOADS=[];
yLOADMODELINFORMATION=[];



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

% CUSTOMER CLASS data
fid2 = fopen([FileName '_CUSTOMERCLASS.txt'], 'w');
tline = fgetl(fid);
while isempty(strfind(tline,'[LOADS]'))
      yCUSTOMERCLASS=tline;
      fprintf(fid2,'%s\r\n',yCUSTOMERCLASS);
      tline = fgetl(fid);
end
fclose(fid2);

% LOADS data
fid2 = fopen([FileName '_LOADS.txt'], 'w');
if ~isempty(strfind(tline,'[LOADS]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_LOADS=SectionID,DeviceNumber,LoadType,Connection,Location'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[CUSTOMER LOADS]')) 
         yLOADS=tline;
         fprintf(fid2,'%s\r\n',yLOADS);
         tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end
 
% CUSTOMER LOADS data
fid2 = fopen([FileName '_CUSTOMERLOADS.txt'], 'w');
if ~isempty(strfind(tline,'[CUSTOMER LOADS]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_CUSTOMERLOADS=SectionID,DeviceNumber,LoadType,CustomerNumber,CustomerType,ConnectionStatus,LockDuringLoadAllocation,Year,LoadModelID,NormalPriority,EmergencyPriority,ValueType,LoadPhase,Value1,Value2,ConnectedKVA,KWH,NumberOfCustomer,CenterTapPercent'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'[LOAD MODEL INFORMATION]')) 
             yCUSTOMERLOADS=tline;
             fprintf(fid2,'%s\r\n',yCUSTOMERLOADS);
             tline = fgetl(fid);
       end
   end    
   fclose(fid2);
end

% LOAD MODEL INFORMATION data
fid2 = fopen([FileName '_LOADMODELINFORMATION.txt'], 'w');
if ~isempty(strfind(tline,'[LOAD MODEL INFORMATION]'))
   tline = fgetl(fid);
   if ~isempty(strfind(tline,'FORMAT_LOADMODELINFORMATION=ID,Name'))
       tline = fgetl(fid);
       while isempty(strfind(tline,'    ')) 
             yLOADMODELINFORMATION=tline;
             fprintf(fid2,'%s\r\n',yLOADMODELINFORMATION);
             tline = fgetl(fid);
             if ~ischar(tline),   break
             end
       end
   end    
   fclose(fid2);
end


ExecutionFlag = 1;

fclose(fid);



    