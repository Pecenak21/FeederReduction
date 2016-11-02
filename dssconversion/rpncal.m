function val = rpncal(expression, varargin)
% evaluate RPN (Reverse Polish Notation) math expressions
% Inputs:
%			expression: expression (string) or cell array of strings to
%			evaluate as rpn commands
% 
% Output:
%			calculated value, or empty value if the input expression was
%			not an RPN value
%
% note: Operators and numbers must be separated by whitespace or commas

if(nargin > 1), expression = [expression varargin]; end
val = [];
s=expression;
while(iscell(s) && length(s)==1)
	s = s{1};
end
	
if ischar(s)
	% clean up expression: remove leading and trailing quotes and spaces
	s = strtrim(s);
	s = regexprep(s,'(^[\(\[\{''"])|([\}\]\)''"]$)','');
	% split on spaces
	s = regexp(s,'[\s,]+','split');
elseif iscellstr(s)
	% remove extra whitespace
	s = strtrim(s);
else
	return;
end


% library of allowed operations
oplib = {'+',...
	'-',...
	'*',...
	'/',...
	'sqrt',...
	'sqr',...
	'^',...
	'sin',...
	'cos',...
	'tan',...
	'asin',...
	'acos',...
	'atan',...
	'atan2',...
	'swap',...
	'rollup',...
	'rolldn',...
	'ln',...
	'pi',...
	'log10',...
	'exp',...
	'inv'};
% check if the given expression is valid RPN
% operation stack
opstack = regexpi(s,'^([a-zA-Z]+[\d]*|[+\-\*\/\^])$','match','once');
isop = ~cellfun(@isempty,opstack);
% for invalid or no operations, return []
if ~any(isop) || any(~ismember(opstack(isop),oplib))
	return;
end
s = str2double(s);
isnum = ~isnan(s);
% if there are any fields that are neither a valid op nor a number, that's
% a fail
if(~all(isop|isnum))
	return;
end
ispi = strcmp(opstack,'pi');
s(ispi) = pi;
isnum = isnum | ispi;

j = 0;
for i=1:length(s)
	% calculate
	if(isnum(i))
		% numbers just get pushed onto the stack
		val = [val; s(i)];
		j = j+1;
	else %isop
		f = str2func(opstack{i});
		switch opstack{i}
			case {'+','-','*','/','^'} % dual input functions
				if(j<2), error('empty RPN stack'); end
				j = j-1;
				val(j) = f(val(j),val(j+1));
			case 'atan2' % dual input, but we have to convert the output from radians to degrees
				if(j<2), error('empty RPN stack'); end
				j = j-1;
				val(j) = f(val(j),val(j+1))*180/pi;
			case {'sqrt','sqr','log10','ln','inv','exp'} % single input functions
				if(j<1), error('empty RPN stack'); end
				val(j) = f(val(j));
			case {'sin','cos','tan','asin','acos','atan'} % single input, but add a 'd' to the function name
				f = str2func([opstack{i} 'd']);
				if(j<1), error('empty RPN stack'); end
				val(j) = f(val(j));
			case 'swap'
				val([j-1 j]) = val([j j-1]);
			case 'rollup'
				val(j+1) = val(j);
				j = j+1;
			case 'rolldn'
				j = j-1;
		end
	end
end

val = val(j);

end

function x = sqr(x)
x = x*x;
end
function x = ln(x)
x = log(x);
end

