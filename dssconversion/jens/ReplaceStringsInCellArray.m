function c = ReplaceStringsInCellArray(c, OldString, NewString)
% c = ReplaceStringsInCellArray(c, OldString, NewString)
%
% PURPOSE : replace string content of a cell with new string
%
%
% INPUT :   c: cell with strings
%           OldString
%           NewString
%
% OUTPUT : cell array with replaced strings

flg=strcmp(c,OldString);
if any(flg)
    c(flg)={NewString};
end

end