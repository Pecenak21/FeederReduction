function [cir cmds] = dssparse(filename)
% Parse OpenDSS file to OpenDSS struct in Matlab
% Outputs:
%			cir : openDSS circuit struct with all components
%			cmds : list of commands/settings to run simulation

%process inputs
if ischar( filename )
	id = find((filename=='/')|(filename=='\'),1,'last');
	fdir = filename(1:id);
	if(~strcmp(filename(1:2),'\\') && filename(2)~=':')
		fdir = [pwd '/' fdir];
		filename = [pwd '/' filename];
	end
else
	error('Invalid input. Must specify filename.')
end

% initialize
cir = struct();
cmds = '';
warningnames = {'dsscapacitor:grounding','cleanPhase:phaselargerthan3'};
for i=1:length(warningnames)
	oldwarnings(i) = warning('off',warningnames{i});
end


% load file
fid = fopen(filename);

ignored = 0;
knownobjs = {};
unknownobjs = {};

while 1
    try
        l = fgetl(fid);
    catch e
        error([filename ' doesn''t exist!']);
    end
    if ~ischar(l),   break,   end
    if strcmp(l,'New Storage.N86 Bus1=SX3104118C.1.2 kV=4.16 kWRated=250 kWhRated=1000 %Reserve=20 kWhStored=1000');
			stopping_ind=1;
	end
	% remove all spaces/tabs at the beginning of the line
	l = strtrim(l);
	
	% skip comments and empty lines; remove any trailing comments
	l = regexp(l,'!|//','split');
	% clean up spaces
	l = regexprep(l,'\s*=\s*','=');
	l = regexprep(l,'\s+',' ');
	
	l = l{1};
	if strcmp(l,'')
		continue;
	end
	% remove 'object=' string if exist
    if ~isempty(strfind(lower(l),'object='))
        i = strfind(lower(l),'object=');
        l(i:i+6) = [];
    end
    
	% Handle the meat content!!!
	cmd = regexp(l,'\S+','match','once');
            
	switch lower(cmd)
		% ignore 'clear' command
		case 'clear'
			continue;
		% handle 'new' command
		case 'new'
            % search for class name and object name
			n  = regexpi(l,'(\S+)\.(\S+)','once','tokens');
			cn = n{1};
            
			if ~ismember(cn,knownobjs) && ~ismember(cn,unknownobjs)
				try
					feval(['dss' lower(cn)]);
					knownobjs = [knownobjs cn];
				catch err
					unknownobjs = [unknownobjs cn];
				end
			end
			
			if ismember(cn,knownobjs)
				ignored = 0;
				obj = createObj(l);
				cn = class(obj);
				cn = cn(4:end);
				if ~isfield(cir,cn)
					cir.(cn) = obj;
				else
					cir.(cn)(end+1) = obj;
				end
			else
				ignored = 1;
			end
		% handle lines start with '~' (continuing of the previous "new" command)
		case '~'
			if ~ignored 
				obj = addtoObj(obj,l);
				cn = class(obj);
				cn = cn(4:end);
				cir.(cn)(end) = obj;
			end
		% handle 'set' command
		case 'set'
			% special handle for base voltages
			if strfind(lower(l),'voltagebases')
				basev = regexp(lower(l),'voltagebases=([\[\(\{"''][^=]+[\]\)\}''"]|[^"''\[\(\{]\S*)','tokens');
				basev = regexp(basev{1}{1},'[\d.]+','match');
				cir.basevoltages = cellfun(@str2num,basev);
			else
				cmds = sprintf('%s\n%s',cmds,l);
			end
		% handle 'redirect' commands
		case {'redirect','compile'}
			fn = regexp(l,'\s+','split');
			fn = strtrim(fn);
			if length(fn) < 2
				error('Check redirect/compile command'); 
			else
				fn = fn{2};
			end
			% check if file name is wrapped in quotes
			m = regexp(l,'[\(\[\{"\''](.*)[\)\]\}"\'']','tokens');
			
			if ~isempty(m)
				fn = m{1}{1};
			end
			
			% get subcir and sub command list
			[cir2 cmdlist] = dssparse([fdir '/' fn]);
			% merge sub circuit to original circuit
			cmds = sprintf('%s\n%s',cmds,cmdlist);
			fnames = fieldnames(cir2);
			for i = 1:length(fnames)
				fn_ = fnames{i};
				if ismember(fnames(i),fieldnames(cir))
					cir.(fn_) = [cir.(fn_) cir2.(fn_)];
				else
					cir.(fn_) = cir2.(fn_);
				end
			end
		% handle other commands
		case 'buscoords' 
			bfn = regexp(l,' ','split');
			bfn = bfn{2};
			fid2 = fopen([fdir '/' bfn]);
			try 
				dat = textscan(fid2,'%[^,]%*[,]%f%*[,]%f%*[\r\n]');
			catch err
				try 
					fseek(fid2,0,'bof');
					dat = textscan(fid2,'%s%f%f%*[\r\n ]');
				catch err
					error('dssparse:buscoordsfileinput','not recoganized input file format. Please use comma or space as delimiter');
				end
			end
			cir.buslist.id = dat{:,1};
			cir.buslist.coord = [dat{:,2}, dat{:,3}];
			fclose(fid2);
		otherwise
			% add to cmds
			cmds = sprintf('%s\n%s\n',cmds,l);
	end
end

for i=1:length(oldwarnings)
	warning(oldwarnings(i).state,oldwarnings(i).identifier);
end

try
	cir.switch = cir.swtcontrol;
	cir = rmfield(cir,'swtcontrol');
catch
end

% close file
fclose(fid);
if ~isempty(unknownobjs)
	disp('Unimplemented object(s):');
	disp(unknownobjs);
end
end

function obj = createObj(l)

% search for class name and object name
[n,prop]  = regexpi(l,'(\S+)\.(\S+)','once','tokens','split');
cn = n{1};
on = n{2};

% create object
obj = feval(['dss' lower(cn)]);
obj.Name = on;

% add properties
if ~isempty(prop{2})
	obj = addtoObj(obj,prop{2});
end

end

function obj = addtoObj(obj,l)

props = regexp(l,'(\S+)=([\[\(\{"''][^=]+[\]\)\}''"]|[^"''\[\(\{]\S*)','tokens');

for i = 1:length(props)
	% clean up quotes
	val = regexprep( props{i}(2),'["'']','' );
	% drop the % notation when given
	prop = regexprep(props{i}(1),'%','');
	% drop the '-' notation when given
	prop = regexprep(prop,'-','');
	
	if strcmp(lower(prop),'bus1')
		Phases=regexp(val,'\.','split');
		Phases=length([Phases{:}])-1;
		if Phases>0
			obj.Phases = Phases;
		end
	end
	
	numval = rpncal(val);
	if isempty(numval)
		obj.(prop) = val;
	else
		obj.(prop) = numval;
	end
end
end