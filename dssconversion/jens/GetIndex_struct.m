function[index] = GetIndex_struct(Y,FieldName,y0)
%function[index] = GetIndex(Structure,FieldName,Content)
%
% PURPOSE : Returns index that contains field content in specified field
%
%
% INPUT : structure, field name, field structure content
%
%
% OUTPUT : index of structure that contains field content
%


index=0;
n=1;

StructureContent=getfield(Y,FieldName);

if isnumeric(y0)
    for i=1:length(StructureContent)
        if StructureContent(i)==y0 
            index(n)=i;
            n=n+1;
        end
    end
else
    for i=1:length(StructureContent)
        if strcmp(StructureContent(i),y0) 
            index(n)=i;
            n=n+1;
        end
    end
end

