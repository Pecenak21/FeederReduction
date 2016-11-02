function[struc_out] = UnifyFormat(struc,FieldName,NewFormat)
% UnifyFormat(struc,FieldName,NewFormat)
%
% PURPOSE : Make the format of the cell ocntent in a structure uniform 
% (i.e., all cells with strings or all cell with numbers)
%
%
% INPUT :   struc: as structure with field names
%           FieldName: Field with the entries that will be unified
%           NewFormat: 'str' or 'num'
%
% OUTPUT :  a structure with reformated entries

n=length(struc.(FieldName));

if strcmp(NewFormat,'str')
    for i=1:n
        if ~ischar(struc.(FieldName){i,1})
            struc.(FieldName)(i,1)={num2str(struc.(FieldName){i,1})};
        end
    end
elseif strcmp(NewFormat,'num')
    for i=1:n
        if ~isnumeric(struc.(FieldName){i,1})
            struc.(FieldName)(i,1)={str2num(struc.(FieldName){i,1})};
        end
    end
end

struc_out=struc.(FieldName);