function s = fnSanitize(s,c)
% Clean up arbitrary strings so that they can be used as struct field names
% Basically we remove anything that's not alphanumeric or _, and prepend an
% 'x' to strings that start with a digit
%
% To activate the prepending 'x' function, must use second argument.
% Example: fnSanitize('042arg','x');
s = strrep(s,'.',' '); %This line is required for UCSD model
s = strtrim(s);
s = regexprep(s,'[^0-9a-zA-Z_]+','_');
s = regexprep(s,'^_+',''); % EMTP-RV needs an alphanumeric first character for bus names
s = regexprep(s,'?','_');

if(nargin >= 2 && ischar(c))
	s = regexprep(s,'^([0-9])',[c '$1']);
end
end



% .          Any single character, including white space
% [c1c2c3]   Any character contained within the brackets: c1 or c2 or c3
% [^c1c2c3]  Any character not contained within the brackets: anything but c1 or c2 or c3
% [c1-c2]    Any character in the range of c1 through c2
% ^expr      Match expr if it occurs at the beginning of the input string.
% expr+      Match expr when it occurs 1 or more times consecutively. Equivalent to {1,}.
 
