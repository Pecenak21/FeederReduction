function x = str2nump(s)
%STR2NUM Convert string matrix to numeric array.
%   X = STR2NUM(S) converts a character array representation of a matrix of
%   numbers to a numeric matrix. For example,
%       
%        S = ['1 2'         str2num(S) => [1 2;3 4]
%             '3 4']
%
%   This works exactly like str2num except that an empty string returns a
%   0 (zero) rather than [].
%   The numbers in the string matrix S should be ASCII character
%   representations of a numeric values.  Each number may contain digits,
%   a decimal point, a leading + or - sign, an 'e' preceding a power of
%   10 scale factor, and an 'i' for a complex unit.
%
%   If the string S does not represent a valid number or matrix,
%   STR2NUM(S) returns the empty matrix.
%
%   CAUTION: STR2NUM uses EVAL to convert the input argument, so side
%   effects can occur if the string contains calls to functions.  Use
%   STR2DOUBLE to avoid such side effects or when S contains a single
%   number.
%
%   Also spaces can be signficant.  For instance, str2num('1+2i') and 
%   str2num('1 + 2i') produce x = 1+2i while str2num('1 +2i') produces
%   x = [1 2i].  These problems are also avoided when you use STR2DOUBLE.
%    
%   See also STR2DOUBLE, NUM2STR, HEX2NUM, CHAR.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.25 $  $Date: 1998/05/12 23:56:38 $

if ~isstr(s) | ndims(s)>2
   error('Requires string or character array input.')
end
%if isempty(s), x = []; return, end
if isempty(s), x = 0; return, end            % this line is by PNath@London.edu
[m,n] = size(s);
if m==1,
  % Replace any char(0) characters with spaces
  s(s==char(0)) = ' ';
  x = protected_conversion(['[' s ']']); % Always add brackets
else
    semi = ';';
    space = ' ';
    if ~any(any(s == '[' | s == ']')), % String does not contain brackets
        o = ones(m-1,1);
        s = [['[';space(o)] s [semi(o) space(o);' ]']]';
    elseif ~any(any(s(1:m-1,:) == semi)), % No ;'s in non-last rows
        s = [s,[semi(ones(m-1,1));space]]';
    else,                               % Put ;'s where appropriate
        spost = space(ones(m,1));
        for i = 1:m-1,
            last = find(fliplr(s(i,:)) ~= space);
            if s(i,n-last(1)+1) ~= semi,
                spost(i) = semi;
            end
        end
        s = [s,spost]';
    end
    x = protected_conversion(s);
end

if isstr(x)
   x = [];
end

function STR2NUM_VaR = protected_conversion(STR2NUM_StR)
% Try to convert the string into a number.  If this fails, return [].
% Protects variables in STR2NUM from "variables" in s.

STR2NUM_LaSTeRR = lasterr;
STR2NUM_VaR = eval(STR2NUM_StR,'[]');
lasterr(STR2NUM_LaSTeRR) % Preseve lasterr