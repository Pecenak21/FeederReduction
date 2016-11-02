function pathtofile = dsswrite( dsscircuit, filename, splitFileFlag, savepath, commands)
% Write OpenDSS circuit to a single file or a set of files with each
% component stored in a file.
% Inputs:
%			dsscircuit: circuit object created from dssconversion function
%			filename: (optional) default: circuit's name. will be used as name for main circuit file and prefix for component files (e.g. [filename]_line.dss )
%			splitFileFlag: (optional) default: 0. Write data to multiple files with each component on a seperate file besides the main one.
%			savepath: (optional) relative/absolute path to save files. If folder doesn't exist, create one.
%			commands: additional commands
% Output:
%			pathtofile: path to main opendss file generated (useful for running OpenDSS Simulation in Matlab)

% Process inputs
if isstruct(dsscircuit)
	c = dsscircuit;
elseif ~isempty(strfind(class(dsscircuit),'dss'))
	splitFileFlag = 0;
	c.(class(dsscircuit)) = dsscircuit;
else
	error('Invalid data type for dsscircuit');
end
% remove wrong fields
cfields = fieldnames(dsscircuit);
for i = 1:numel(cfields)
	if isempty(strfind(class(dsscircuit.(cfields{i})), 'dss')) && ~strcmpi(cfields{i},'buslist') && ~strcmpi(cfields{i},'basevoltages')
		c= rmfield(c,cfields{i});
	end
end

if ~exist('splitFileFlag','var')
	splitFileFlag = 1;
end

headerfooterflag = 1;

if ~exist('savepath','var') || isempty(savepath)
	savepath = [pwd '/tmp/'];
    if ~exist(savepath,'dir'), mkdir(savepath); end
else
	% handle relative path
	if(~strcmp(savepath(1:2),'\\') && savepath(2)~=':')
		savepath = [pwd '/' savepath];
	end
	
	% create folder if it doesn't exist
	if exist(savepath,'dir') < 1
		mkdir(savepath);
	end
end

% handle filename. Use circuit name if not specified
if exist('filename','var') && ~isempty(filename)
	if strfind(filename,'.dss'), 
		fname = filename(1:strfind(filename,'.dss')-1); 
	else
		fname = filename;
	end
else
	if ~isfield(c,'circuit') 
		warning('dsswrite:circuitUndefined','The input data doesn''t contain a circuit object! Check and make sure this is what you want.\nOpenDSS will not be able to load this file by itself.');
		fname = 'newdssfile';
		headerfooterflag = 0;
	else
		fname = c.circuit.Name;
	end
end

if exist([savepath '/' fname '.dss'],'file')
	try
		fnn = [savepath '/' fname '.dss'];
		fid = fopen(fnn,'w');
		fclose(fid);
	catch
		fname = [fname '_' datestr(now,'YYYYMMDDhhmmss')];
	end
end

% Writing out
s = '';
if isfield(c,'circuit')
	s = char(c.circuit);
end
fn = fieldnames(c);
% bus list isn't a class like the others
ind = strcmp(fn,'buslist')|strcmp(fn,'basevoltages');
if(any(ind))
	fn(ind) = [];
	if isfield(c,'buslist')
		buslist = c.buslist;
	end
end
% arrange them in the right order to print out (e.g. linecode should be
% defined before line)
classes = {'wiredata','linegeometry','linecode','line','loadshape','tshape','tcc_curve','reactor','fuse','transformer','regcontrol','capacitor','capcontrol','xycurve','pvsystem','InvControl','storage','storagecontroller','load','generator','swtcontrol','monitor','energymeter'};
classes(~ismember(classes,fn)) = [];

classes2 = {'switch','recloser'};
classes2(~ismember(classes2,fn)) = [];

fn = [classes, setdiff(fn,[classes classes2])', classes2];

% open main file
if isempty(strfind(fname,'.dss')), fname2 = [fname '.dss']; end
pathtofile = [savepath '/' fname2];
fidmain = fopen(pathtofile, 'w');
if fidmain < 1
	splitFileFlag = 0;
	fname2 = ['dss_' datestr(now,'YYYYMMDDhhmmss')];
	pathtofile = [savepath '/' fname2 '.dss'];
	fidmain = fopen(pathtofile, 'w');
end

if ~splitFileFlag
	% add each device class to circuit string
	for i = 1:length(fn)
		if strcmp('circuit',fn{i}), continue; end;
        if ~isempty(c.(fn{i}))
            s = [s char(c.(fn{i}))];
        end
	end
else
	% write files for all devices
	
	for i = 1:length(fn)
		
		%device filename
		dfn = [fname '_' fn{i} '.dss'];
		
		if strcmp('circuit',fn{i}), continue; end;
		% open file for writing
		fid = fopen([savepath '/' dfn], 'w');
		if(fid==-1), error('dsswrite:openfailed','Failed to open output file %s for writing!\nRemember to close open files before overwriting.',[savepath '/' dfn]); end
		s_ = char(c.(fn{i}));
		try
			fwrite(fid, s_);
			fclose(fid);
		catch err
			warning('dsswrite:openfiles','Remember to close files before overwriting them!');
			rethrow(err);
		end
			
		% update main file with "Redirect" command
		s = sprintf('%s\n%s',s,['Redirect ' dfn]); 
	end
end

if(~isfield(c,'basevoltages'))
	c.basevoltages = [115, 12.47, 4.16, 0.48, 0.12];
end

if headerfooterflag
	cvs = sprintf('\n\n! Let DSS estimate the voltage bases\n%s%s\n%s\n',...
		'Set voltagebases=',mat2str(c.basevoltages),...
		'Calcvoltagebases     ! This also establishes the bus list');
if isfield(c,'InvControl')
	iter=sprintf('\nSet maxcontroliter=5000');
	s = [s cvs iter];
else
	s = [s cvs];
end
	
	
end


	

if(exist('buslist','var'))
	s = sprintf('%s\nBuscoords %s_%s.csv\n', s, fname, 'buscoords');
	sbl = [buslist.id num2cell(buslist.coord)]';
	sbl = sprintf('%s, %g, %g\n',sbl{:});
	try
		fid = fopen([savepath '/' fname '_buscoords.csv'], 'w');
		fwrite(fid,sbl);
		fclose(fid);
	catch err
		warning('dsswrite:openfiles','Remember to close files before overwriting them!');
		rethrow(err);
	end
end

if headerfooterflag
	% header
	h = sprintf('%s\n\n','Clear');
	
	% footer
% 	f = sprintf('\n\n%s\n%s\n\n%s\n%s\n%s\n',...
% 				'set maxiterations=100',...
% 				'solve mode=snapshot',...
% 				'show voltages LL Nodes',...
% 				'show powers kva elements',...
% 				'show taps',...
%                 'export voltages',...
%                 'export seqvoltages',... 
%                 'export powers kva',...
%                 'export p_byphase',...
%                 'export seqpowers');		
	f = sprintf('\n\n%s\n%s\n\n%s\n%s\n%s\n',...
				'set maxiterations=1000');
else
	h = '';
	f = '';
end

% write main file
if ~exist('commands','var')
	fwrite(fidmain, [h s f]);
else
	fwrite(fidmain, sprintf('%s %s \n %s',h,s,commands));
end
fclose(fidmain);

end