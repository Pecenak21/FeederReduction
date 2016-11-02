function [d, newname, rawd] = importRedlandsData(file)
% d = importRedlandsData can be used for reading the very large xlsx files that we recieve from SCE
% 
% recommended usage is something like this
%
% idir = '/path/to/idir';
% odir = '/mnt/lab_18tb1/database/deployments/Redlands/matlab_data';
% flist = dir([idir '/*.xlsx']);
% for i=1:length(flist)
%	[d, n] = importRedlandsData([idir '/' flist(i).name]);
%	save([odir '/' n], 'd');
% end
%
% Note: this function requires gnu sed for large files, which probably means you can't run it on windows!

% basically limited reading of an xlsx file.  A lot of things are stolen from Matlab's built in routines
% I'm cobbling this together so as to avoid needing to wait for matlab to swap itself into oblivion when it tries to parse the whole file at once

% calculate the new filename
[~,newname] = fileparts(file);
try
	% extract the parts
	nn = regexp(newname,'([^. -]*)','tokens');
	nn = horzcat(nn{:});
	dmask = ~cellfun(@isempty,regexp(nn,'\d{6}','once'));
	nn = horzcat(nn(dmask),nn(~dmask));
	nn(2,:) = {'_'}; nn{2,end} = '.mat';
	newname = lower(strcat(nn{:}));
catch %#ok<CTCH>
	newname = [lower(regexprep(newname,'[ -]+','_')) '.mat'];
end

% for small files, we'll use the built-in xlsread.  This has the advantage of working with non-rectangular data files, but the disadvantage of not working with very large files
% I implemented it this way because of running into non-rectangular tables in some small files, but if you hit any of them in large files, you're probably just SOL, or you'll have to figure out how to make all this work the other way...
f_inf = dir(file);
if f_inf.bytes <= 20e6 % file size less than 20MB
	fprintf('Small file; calling xslread()...'); t = tic;
	[~,~,d_] = xlsread(file);
	toc(t); fprintf('Extracting data...'); t = tic;
	% get the row header data
	header = d_(1:2,1:4:end);
	d_(1:2,:) = [];
	
	% extract the data; comes in groups of four columns
	for i=1:size(header,2)
		% description is the long line of the header
		d(i).desc = header{1,i};
		% index of first data column:
		ii = i*4-3;
		% index of first empty row for this field:
		m = find(~cellfun(@ischar,d_(:,ii)),1,'first');
		if(isempty(m)) % this column goes all the way to the end
			m = size(d_,1)+1;
		end
		% timestamps + values
		d(i).time = datenum(d_(1:m-1,ii)) + vertcat(d_{1:m-1,ii+1});
		d(i).v = vertcat(d_{1:m-1,ii+3});
		% quality control flags
		qcf = regexprep(d_(1:m-1,ii+2),' ','');
		d(i).qc = uint8(~strcmpi(qcf,'ok'));
		m = ~d(i).qc;
		mm = m & cellfun(@isempty,regexp(qcf,'^[\.#]*$'));
		d(i).qc(mm) = cellfun(@(x)sum(x(:)=='#'),qcf(mm))+1;
		if(any(m & ~mm))
			mm = m&~mm;
			% this accomplishes more or less the same thing as grouping unique used below, but works for strings and not quite as fast.  We could also just bite the bullet and re-modify a copy of unique to make a more-global grouping_unique.
			qcf_u = unique(qcf(mm));
			[~,j] = ismember(qcf(mm),qcf_u);
			d(i).qc(mm) = numel(qcf_u)+j;
		end
		
		for j = 1:length(qcf)
			if(strcmpi(qcf{j},'ok'))
				qcf{j} = 0;
			elseif(regexp(qcf{j},'^[\.#]*$'))
				qcf{j} = sum(qcf{j}=='#')+1;
			else
				qcf{j} = length(qcf)+j;
			end
		end

	end
	toc(t);
	
	% skip the painful method of reading the files
	return;
end

% unzip, load the strings, and get a file path for the sheet we want, essentially stolen verbatim from Matlab
baseDir = tempname;
cleanupBaseDir = onCleanup(@()rmdir(baseDir,'s'));
fprintf('Unzipping %s...',file); t = tic;
unzip(file, baseDir);
toc(t);

sharedStrings = extractSharedStrings(baseDir);

sheetIndex = 1;
workSheetFile = fullfile(baseDir, 'xl', 'worksheets', sprintf('sheet%d.xml', sheetIndex));
workSheetTSV = regexprep(workSheetFile,'xml$','tsv');

% Run sed to pare it down to just the data:
fprintf('Extracting TSV data...'); t = tic;
% use matlab to put in newlines for the rows
% =======================
fidx = fopen(workSheetFile);
rawd = fread(fidx,'*char')';
fclose(fidx);

ridx = strfind(rawd,'<row');
nl = sprintf('\n');
rawd(ridx) = nl;
fidx = fopen([workSheetFile '_lines'],'w');
fwrite(fidx,rawd);
fclose(fidx);
% =======================
unix(['sed -r -e ''1,2d'' -e ''s/^row[^>]*>//g'' -e ''s/<\/?[c][^>]*>//g'' -e ''s/^<v>//'' -e ''s/<v>/	/g'' -e ''s/<\/v>//g'' -e ''s/<\/row.*//'' <' workSheetFile '_lines >' workSheetTSV]);
toc(t);
[~, data_len] = unix(['wc -l ' workSheetTSV]);
data_len = sscanf(data_len,'%d');

% open the file for reading
fprintf('Reading and formatting data...'); t = tic;
fid = fopen(workSheetTSV);
cleanupFH = onCleanup(@()fclose(fid));

% The first two lines of the resulting TSV are headers.  Read them:
header = regexp([fgetl(fid) '	' fgetl(fid)],'\t','split');

% read the rest of the data as numeric types
n = ftell(fid);
x = sscanf(fgetl(fid),'%f')';
fseek(fid,n,'bof');
rawd = fscanf(fid,'%f',[numel(x),inf])';
if(data_len-size(rawd,1)~=1) % two header rows, minus 1 for wc messing up because there's not a trailing newline means the difference should be 1
	warning('importRedlands:headerNoMatch','\nfile does not appear to contain rectangular data.  Expected %d rows of data, found %u.\n', data_len-1,size(rawd,1));
end

% straighten up the header
header = str2double(header);
header = reshape(header,numel(header)/2,2)';
header = sharedStrings(header+1);

% data comes in 4-column blocks
%	column 1 is the day
%	column 2 is the fractional day
%	column 3 seems to be a QC field
%	column 4 is the useful data
% each block gets two rows worth of header

d = struct('time',cell(1,size(header,2)), 'qc',[], 'v',[]);
for iG=1:size(header,2)
	i = iG*4-3;
	if(i>size(rawd,2))
		fid2 = fopen(workSheetFile);
		x = fread(fid2,10000,'*char');
		fclose(fid2);
		x = regexp(x','<row','split');
		x = regexp(x,'<c\s+r="(?<ranges>[A-Z]+\d+)"','names')';
		headercols = regexprep({x{2}.ranges},'1$','');
		datacols = regexprep({x{4}(1:4:end).ranges},'3$','');

		warning('importRedlands:headerNoMatch','The header in this file does not have the same number of columns as the data fields.  I may not have matched them correctly.');
		break;
	end % can't fill in the last several columns
	% fill in dates in some of the columns
	[dates,ia] = grouping_unique(rawd(:,i));
	dates = datenum(sharedStrings(dates+1));
	for j=1:length(ia)
		rawd(ia{j},i) = dates(j);
	end
	% fill in qc flags in column 3
	[qcf,ia] = grouping_unique(rawd(:,i+2));
	qcf = regexprep(sharedStrings(qcf+1),' ','');
	for j = 1:length(qcf)
		if(strcmpi(qcf{j},'ok'))
			qcf{j} = 0;
		elseif(regexp(qcf{j},'^[\.#]*$'))
			qcf{j} = sum(qcf{j}=='#')+1;
		else
			qcf{j} = length(qcf)+j;
		end
	end
	qcf = cell2mat(qcf);
	d(iG).qc = uint8(size(rawd(:,i+2)));
	for j=1:length(ia)
		d(iG).qc(ia{j}) = qcf(j);
	end

	% copy over time and data fields
	d(iG).time = rawd(:,i)+rawd(:,i+1);
	d(iG).v = rawd(:,i+3);
end
% realign data to go with descriptions
if(exist('datacols','var'))
	[ix,x] = ismember(datacols,headercols);
	d(x) = d(ix);
	ix = ~ismember(headercols,datacols);
	d(ix) = struct('time',[], 'qc',[], 'v',[]);
end
% assign descriptions
[d.desc] = deal(header{1,:});

toc(t);

end

% This one is stolen verbatim from matlab in its entirety...
function sharedStrings = extractSharedStrings(baseDir)

fid  = fopen(fullfile(baseDir, 'xl', 'sharedStrings.xml'), 'r', 'n', 'UTF-8');
sharedStrings = '';
if fid ~= -1
	sharedStrings = fread(fid, 'char=>char')';
	fclose(fid);
	
	% Rich text is captured across multiple "<t>" nodes (so the number of <t> nodes may outnumber
	% the of number of <si> (the real number of strings), we need to concatinate mulitiple <t>'s
	% into single <si>'s.
	stringItemElements = regexp(sharedStrings,'<si>(.*?)</si>','tokens');
	
	sharedStrings = {};
	if ~isempty(stringItemElements)
		groupedTextElements = regexp([stringItemElements{:}],'<t.*?>(?<textElements>.*?)</t>','names');
		for i=length(groupedTextElements):-1:1
			if isempty(groupedTextElements{i})
				sharedStrings{i} = '';
			else
				sharedStrings{i} = [groupedTextElements{i}.textElements];
			end
		end
	end
end
end
