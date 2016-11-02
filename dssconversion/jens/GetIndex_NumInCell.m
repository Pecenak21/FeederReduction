function[index, y] = GetIndex_NumInCell(cell,num)
% function[index y] = GetIndex_StrInCell(cell,str)
%
% PURPOSE : returns indeces of str in cell array of strings, exact match
%
%
% INPUT :  cell: cell with numbers
%           num: number to find
%
% OUTPUT : index: array with locations of num
%              y: array with 1s at locations that matches num


num_array=cell2mat(cell);
index=find(num_array==num);
if ~isempty(index)
    y=index(1);
else
    y=[];
end
