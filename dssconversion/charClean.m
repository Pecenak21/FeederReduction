%% handle char input
function val = charClean(val)
% output: either a single simple string or an clean cellstring

if iscell(val) && length(val) == 1
	val = val{1};
end

% clean up string
% trim space at the beginning and end
val = strtrim(val);

% ASSUMPTION: if string represent a matrix, rows will be seperated by | or semicolon (not comma). Values in a row are seperated by either comma or space.
% clean up unnecessary spaces and replace comma by single space
val = regexprep(val,'[\s\,]+',' ');
% clean up spaces around semicolon, | if represented and replace by | |
val = regexprep(val,'\s*[\|;]\s*','\|');
% clean up spaces around matrix brackets
val = regexprep(val,'\s*[\]\)\}]\s*','\}');
val = regexprep(val,'\s*[\[\(\{]\s*','\{');
% clean up quotes
val = regexprep(val,'"''','');

% three possibilities are considered here: matrix of num, cell of
% strings and normal string

% ASSUMPTION: matrix is surrounded by [] or {} or () and array of
% strings will be 
%-> This Regexp has a bug in it.  E.g. try val = '[abc 123] [xyz]'.
%-> Not fixing it now because I don't think it's important
m = regexp(val,'[\(\[\{](.*)[\)\]\}]','tokens');
if ~isempty(m)
	% it does represent a matrix
	if length(m) > 1 % contain multiple matrices
		warning('dataclean:multiMatrixStringInput','String contain multiple matrices. Only process the first matrix. Ignore the remaining!');
		val = m{1}{1};
	else %length(m) == 1: single matrix (yahhooo)
		if strcmp(m,' ')
			return;
		elseif ~any(val=='|')
			val = regexp(m{1}{1},'\s','split');
		else
			val = m{1}{1};
		end
	end
end
% seperate string into cellstring if it contains space or |
if ~iscell(val) && length( regexp(val,'[\|;]','split') ) > 1
	% split by row
	x = regexp(val,'[\|;]','split')';
	% split by colume
	for i = 1:length(x)
		if isempty(x{i}), continue; end;
		y = regexp(x{i},'[\,\s]','split');
		for j = 1:length(y)
			newval{i,j} = y{j};
		end
	end
	val = newval;
end

if iscell(val) && length(val) == 1
	val = val{1};
end
		
end