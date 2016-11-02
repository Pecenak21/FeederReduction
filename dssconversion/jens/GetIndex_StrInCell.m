function[index, y] = GetIndex_StrInCell(cell,str,k)
% function[index y] = GetIndex_StrInCell(cell,str)
%
% PURPOSE : returns indeces of str in cell array of strings, exact match
%
%
% INPUT :  cell: cell with strings
%           str: string to find
%             k: first k hits are returned
%
% OUTPUT : index: array with string locations
%              y: array with 1s at locations that matches str

if nargin<3
    k=0;
end

y=strcmp(cell,str); %alternatively, the strfind function can be used

if k>0
    index=find(y,k);
else
    index=find(y);
end

% this is a brute force find, there must be a smarter way to determine a
% non-zero cell in a cell array
% perhaps that: matdata=cellfun(@str2num,data),
% find(~cellfun(@isempty,ans))
n=1;
if isempty(index)
    y=strfind(cell,str);
    for i=1:length(y)
        if ~isempty(y{i})
            if n==1 % only returns and index if unique hit
                index=i;
                n=n+1;
            else
                index=[];
            end
        end
    end
end
