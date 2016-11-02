function[y index] = FindStrInCell(cell,str,flgExact)
% function[y] = FindStrInCell(cell,str)
%
% PURPOSE : returns 1 if string is a substring of strings in a cell
%
%
% INPUT :  cell: cell with words
%          str: substring
%          flgExact: exact match required, default is 0
%
% OUTPUT : y: 1 if str is a substring in the cell strings
%             0 if str is not a substring in any of the cell strings

if nargin<3
    flgExact=0;
end

if isempty(cell)
    y=0;
    return
end

for i=1:length(cell);
    if flgExact
        if strcmp(cell{i},str)
            y=1;
            index=i;
            return
        else
            y=0;
            index=0;
        end
    else
        if ~isempty(strfind(cell{i},str))
            y=1;
            index=i;
            return
        else
            y=0;
            index=0;
        end
    end
end






