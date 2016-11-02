function [v varargout]= dataclean(val,type,outputtype)
% [v varargout] = dataclean(val,type,outputtype)
% clean up data. An effort to capture as many as possible input
% types/representations. Assumptions are made based on most common
% conventions to present data. 
%
% Input:
%			val: input value
%			type: intended type of input value. 
%				Supported type: 'num','logical','name','phase','conn','monitoredphase'
%			outputtype: (optional) default: same as 'type'. If specified, output will be convert to this type.
% Output:
%			v: output value with type 'outputype' if specified. Otherwise will have type 'type'
%			varargout: additional output for special cases
%				[conn, grounded] = dataclean(val,'conn')
%				[phasenum phasestring phaseOpenDSSform] = dataclean(val,'monitoredphase')
% 
% Examples of use:
%			v = dataclean('[1 | 2 3 | 4 5 6]','num')
%			v = dataclean('TrUe','logical') % output 1
%			v = dataclean('TrUe','logical','char') % output 'TRUE' instead of 1
%			[conn, grounded] = dataclean({'wye','del'},'conn')
%			[phasenum phasestring phaseOpenDSSform] = dataclean('[xyzg xzg]','monitoredphase')
%			
%
% Assumptions:
%		1. if string represent a matrix, rows will be seperated by | or semicolon (not comma). Values in a row are seperated by either comma or space.
%		2. matrix is surrounded by [] or {} or (). If not, it should be		seperated by | or ; and spaces
% 
% NOTES: outputtype option only works for 'logical' type only now.

% initialize
v = []; 

% process inputs
if ~exist('val','var') || isempty(val)
	return;
end
if ~exist('type','var')
	type = '';
end
if ~exist('outputtype','var')
	outputtype = '';
end

%% determine data type from input value
% find out what kind of data we're having here
% cleanup if it's a multi-layer cell
while iscell(val) && length(val) == 1
	val = val{1};
end

% preprocess data
% now val can only be either a char, a cellstring, a cell array, a struct, a function handle, a logical or a numeric value
% check if val is an array of value in string format
% e.g. [1 2 3]
if iscell(val)
	% cell string
	if ischar(val{1})
		if length(val) > 1
			val = cellfun(@charClean,val,'UniformOutput',0);
		else
			val = charClean(val{1});
		end
	elseif isnumeric(val(1)) || islogical(val(1))
		% cell array
	else
		% cell of structs, function_handle(s) or custom class
	end
elseif ischar(val)
	val = charClean(val);
elseif isnumeric(val)
elseif isa(val,'struct')
elseif isa(val,'function_handle')
else
	% custom class
end

% handle cell string inputs
iscs = iscellstr(val);

switch lower(type)
	case 'name'
		v = cleanName(val);
	case 'num'
		if iscs
			v = cellfun(@cleanNum,val);
		else
			v = cleanNum(val);
		end
	case 'phase'
		if iscs
			v = cellfun(@cleanPhase,val);
		else
			v = cleanPhase(val);
		end
	case {'conn','connection'}
		if iscs
			[v varargout{1}] = cellfun(@cleanConn,val,'UniformOutput',0);
		else
			[v varargout{1}] = cleanConn(val);
		end
	case {'monitoredphase','monphase'}
		if iscs
			[v varargout{1} varargout{2}] = cellfun(@cleanMonitoredPhase,val,'UniformOutput',0);
		else
			[v varargout{1} varargout{2}] = cleanMonitoredPhase(val);
		end
	case 'logical'
		if iscs
			v = cellfun(@cleanLogical,val,outputtype);
		else
			v = cleanLogical(val,outputtype);
		end
	case 'file'
		warning('cleandata:tobeimplemented', 'Return input value for now.');
		v = val;
	otherwise
		warning('cleandata:datatype','not supported data type. Return input value.');
		v = val;
end

end

%% cleanName
function val = cleanName(val)
val = regexprep(val,'[ ''"-]+','_');
end

%% cleanLogical value
function val = cleanLogical(val,outputtype)

try 
	if islogical(val)
	elseif isnumeric(val)
		val = logical(val);
	elseif ischar(val)
		switch lower(val(1))
			case {'0','n','f'}
				val = false;
			case {'1','y','t'}
				val = true;
			otherwise
				warning('dataclean:cleanLogical','not supported input for logical type. Set to 0.');
				val = 0;
		end
	end
catch err
	warning('dataclean:cleanLogical','not supported input for logical type. Set to 0.');
	val = false;
end

% output as char when requested
if exist('outputtype','var') && ismember(lower(outputtype),{'string','char'})
	if val
		val = 'TRUE';
	else
		val = 'FALSE';
	end
end

end
%% cleanNum
function val = cleanNum(s)
% clean up numeric inputs. Handle double, string (containing single value
% or array of number in any quotes), and cellstring with each cell is a number.

% initialize output
val =[];

if ~exist('s','var')
	warning('cleandata:cleanNum','You must provide an input');
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
elseif iscell(s)
	val = str2double(s);
	if any(isnan(val(:)))
		val = genfullmatrix(val);
	end
elseif islogical(s)
	val = s;
else
	warning('cleandata:cleanNum','Invalid input for numeric values!');
	val = s;
end
end

%% cleanConn
function [conn, grounded] = cleanConn(str)
% clean up connection for OpenDSS

% clean input
str = regexprep(str,'[\s"'']','');

% check if input is in array format
id = regexp(str,'[\[\(\{\]\)\}]');
if ~isempty(id)
	val(id) = '';
	val = regexp(val,'[,\s]','split');
	for i = 1:length(val)
		val{i} = getconn(val{i});
	end
else
	conn = getconn(str);
end

% Check if grounded
if strcmp(conn,'wye') % only consider grounded when wye conn is used.
	grounded = double(any(lower(str)=='g'));
else
	grounded = 0;
end

end

function conn = getconn(str)
switch lower(str(1))
	case {'w','y'} % {'wye','yg','y'}
		conn = 'wye';
	case {'d'} % {'delta','d','del'}
		conn = 'delta';
	case {'l'} % two cases here: 'ln' (wye), 'll' (delta)
		if strcmp( str(2),'n')
			conn = 'wye';
		else 
			conn = 'delta';
		end
	otherwise
		warning('dataclean:getconn','invalid connection type. Set to wye connection.');
		conn = 'wye'; return;
end
end

%% cleanPhase
function [ numPhase ] = cleanPhase( phase )
%PHASECLEAN Clean up phase input and return number of phases for OpenDSS
%conversion. Default numPhase: 3.

% in case s is represented as a string, e.g '3'
try 
	s_ = str2double(phase);
	if(1<=s_)
		numPhase = s_;
		if(s_>3)
			warning('cleanPhase:phaselargerthan3','There are more than 3 phases.');
		end
		return;
	end
catch err
	% it's not, keep going
end

% process input
if isempty(phase)
	warning('DSSConversion:phaseClean','Need to specify input. Set to 3.');
	phase = 3;
end

if ischar(phase)
	% remove spaces, neutral and ground notation
	phase = regexprep(lower(phase),'[ng\s]','');
	
	% count number of phases excluding neutral line
	phase = length(phase);
elseif ~isnumeric(phase)
	warning('DSSConversion:phaseClean','Invalid input for phases. Set to 3.');
	phase = 3;
end

numPhase = floor(phase);

if numPhase < 0
	warning('DSSConversion:phaseClean','Invalid numeric input for phases. Set to 3.');
	numPhase = 3;
end
if numPhase > 3
	warning('DSSConversion:phaseClean',['NumPhase = ' num2str(numPhase) ' > 3: More than 3 phases.']);
end

end

%% cleanMonitoredPhase
function [p ps pfullform] = cleanMonitoredPhase(s)
% clean up mornitored phase.
% output:
%			p : {1,2,3} phase number; for delta/ll connection use the first of the two phases (so 1 for 1-2, 2 for 2-3, 3 for 3-1). Default to 1 if input is invalid.
%			ps: phase(s) monitored in form of 'xyz'
%			pfullform: mornitored phase represented as full form for OpenDSS bus. E.g: 1.2.3, 1.2, 2.3

% in case s is represented as a string, e.g '3'
try 
	s_ = str2double(s);
	if(1<=s_ && s_<=3)
		s = s_;
	end
catch err
	% it's not, keep going
end

if ischar(s)
	s(s==' ') = '';
	s = sort(lower(s));
	allphase = 'abcxyz';
	[val id] = ismember(s,allphase);
	ps = allphase(nonzeros(id));
	ps(ps=='a') = 'x';
	ps(ps=='b') = 'y';
	ps(ps=='c') = 'z';
	ps = unique(ps);
	switch ps
		case {'x','xy'}
			p = 1;
		case {'y','yz'}
			p = 2;
		case {'z','xz'}
			p = 3;
		case {'xyz'} % Using all three phases (p should be ignored)
			p = 0;
		otherwise
			warning('DSSConversion:mornitoredPhaseClean','Invalid monitored phase input');
			p = 1; ps ='x';
	end
elseif isnumeric(s)
	p = uint8(s);
	switch p
		case 1
			ps = 'x';
		case 2
			ps = 'y';
		case 3 
			ps = 'z';
		otherwise
			warning('DSSConversion:mornitoredPhaseClean','Invalid monitored phase input');
			p = 1; ps = 'x';
	end
else
	warning('DSSConversion:mornitoredPhaseClean','Invalid monitored phase input');
	p = 1; ps = 'x';
end
	pfullform = generateFullForm(ps);
end

%% general functions
function o = generateFullForm(ps)
	if length(ps) == 1, o = num2str(toNum(ps));
	elseif length(ps) > 1
		o = num2str(toNum(ps(1)));
		for i=2:length(ps)
			o = [o '.' num2str(toNum(ps(i)))];
		end
	else
		error('Invalid input');
	end
		
end

function o = toNum(p)
	switch p
		case {'a','x'}
			o = 1;
		case {'b','y'}
			o = 2;
		case {'c','z'}
			o = 3;
	end
end

function v = genfullmatrix(d)
% works for lower matrix only now
if ischar(d)
	% handle string in form {'1' '2 3' '4 5 6'}
	l = length(d);
	v = nan(l);
	for i = 1:l
		v(i,1:i) = d{i};
	end
elseif isnumeric(d)
	% handle numeric lower triangle matrix input
	v = d;
end

% fill the upper half
for i = 1:size(d,1)
	for j = 1:size(d,1)
		if isnan(v(i,j)), v(i,j) = v(j,i); end;
	end
end
	
end