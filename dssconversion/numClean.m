function val = numClean(s)
% clean up numeric inputs. Handle double, string (containing single value
% or array of number in any quotes), and cellstring with each cell is a number.

% initialize output
val =[];

if ~exist('s','var')
	warning('You must provide an input');
	return;
end

if isnumeric(s)
	val = s;
elseif ischar(s)
	qid = regexp(s,'[\[\]\{\}\(\)"'''']');
	if ~isempty(qid)
		s(qid)=''; 
	end
	% handle matrix data
	% can only handle form: [1 2 3| 4 5 6] at the moment
	if ~isempty(regexp(s,'\|','once'))
		num = regexp(s,'\|','split');
		val = cell(1,length(num));
		for i = 1:length(num)
			val{i} = sscanf(num{i},'%f')';
		end
		if length(val{1})~=length(val{2})
			val = genfullmatrix(val);
		end
	else
		val = sscanf(s,'%f')';
		if isempty(val)
			val = s;
		end
	end
elseif iscellstr(s)
	val = str2double(s);
	if any(isnan(val))
		val = s;
	end
elseif islogical(s)
	val = s;
else
	warning('Invalid input for numeric values!');
	val = s;
end
end

function v = genfullmatrix(d)
% handle [1|2 3| 4 5 6]
% works for lower matrix only now
l = length(d);
v = nan(l);
for i = 1:l
	v(i,1:i) = d{i};
end
% fill the upper half
for i = 1:l
	for j = 1:l
		if isnan(v(i,j)), v(i,j) = v(j,i); end;
	end
end
end