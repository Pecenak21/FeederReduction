function obs = cymeread(fns)
% ob = cymeread(fns) reads CYME formatted data from the files listed in cell array fns
%	cymeread has the following behaviors:
%		* each input file becomes a field of the return value
%		* each section of each input file becomes a field of that file's struct
%		* sections with a single header line are represented as a struct array with fieldnames corresponding to the header line (aka "table struct")
%		* sections with two named formats become a struct with two fields which are each a "table struct" corresponding to their parts of the file
%		* sections with two different formats of the same name become a cell array with "table struct"s containing the data
	for i_ = 1:length(fns)
        clear ob;
		fn=fns{i_};
		if(~exist(fn,'file')), continue; end
		
		% dump the whole file into memory for dealing with conveniently
		flines = regexp(fileread(fn),'\r?\n','split');
		
		% breaks between sections are blank lines, and sections start on the next line after that
		sectionbreaks_2 = ~cellfun(@isempty,regexp(flines,'^\s*$','once'));
		sectionbreaks = find(cellfun(@isempty,flines) | sectionbreaks_2);
		sectionstarts = [1 sectionbreaks+1];
		mask = [diff(sectionstarts)==1 sectionstarts(end)>numel(flines)];
		sectionstarts(mask) = [];
		
		% remove the blank lines so we don't have to worry about them.
		for sn = 1:length(sectionstarts)
			sectionstarts(sn) = sectionstarts(sn) - sum(sectionbreaks<sectionstarts(sn));
		end
		flines(sectionbreaks) = [];
		
		% loop over the sections and process them
		sectionends = [sectionstarts(2:end)-1 numel(flines)];
		for sn = 1:length(sectionstarts)
			sec_name = regexp(flines{sectionstarts(sn)},'^\[(.*)]','once','tokens');
			if(isempty(sec_name)); warning('cymeread:corruptfile','sections must start with a [] block, but the first line of this one is ''%s''.  Skipping it.',flines{sectionstarts(sn)}); continue; end
			sec_name = sanitize(lower(sec_name{1}));
			if(strcmp(sec_name,'general'))
				ob.(sec_name) = dogen(flines(sectionstarts(sn)+1:sectionends(sn)));
			else
				ob.(sec_name) = dodata(flines(sectionstarts(sn)+1:sectionends(sn)));
			end
		end
		
		[fn fn] = fileparts(fn);
		fn = sanitize(lower(regexpi(fn,'equipment|load|network','once','match')));
		obs.(fn) = ob;
	end
	
end

function g = dogen( lines )
% Read in the GENERAL section, which is different because it has
% lines like:
%	DATE=MAY 16,2012
%	CYME_VERSION=xxx
f = regexp(lines,'=','once','split');
f = vertcat(f{:})';
f(1,:) = lower(f(1,:));
g = struct(f{:});
end

function s = dodata( lines )
% Read in other sections, which have a header line
%	FORMAT_XXX=fieldname1,fieldname2,...
% followed by comma-delimited data lines
s = [];
% try to read the header line
% the case where it's missing is generally the unit specification (e.g. IMPERIAL)
if(isempty(lines)), return; end;

cData = regexp(lines,',','split');
%firstF = cellfun(@(x)x{1},cData,'UniformOutput',0);
specialF = regexp(lines,'^(\w*)=','once','tokens'); % ~cellfun(@isempty,regexp(firstF,'='));
specialI = ~cellfun(@isempty,specialF);
formatI = ~cellfun(@isempty,regexp(lines,'^FORMAT'));
if(sum(formatI)==1)
	% first line is format; no others should be special.  this is the easy case
	if(sum(specialI)>1); warning('cymeread:unknownspecial','Unknown special line in this section'); end
	dat = vertcat(cData{2:end});
	sfns = cData{1};
	sfns{1} = regexprep(sfns{1},'^\w*=','');
	sfns = sanitize(sfns);
	s = cell2struct(dat,sfns,2);
elseif(all(find(formatI)<=sum(formatI)))
	% first n lines are formats, subsequent special lines use one of the other formats
	specials = lower(regexprep(cellfun(@(x)x{1},cData(specialI),'UniformOutput',0),'^(\w*)=.*','$1'));
	formatI = find(formatI);
	gfns = regexprep(cellfun(@(x)x{1},cData(formatI),'UniformOutput',0),'^FORMAT_(\w+)=.*','$1');
	gfns = sanitize(lower(gfns));
	for i=1:length(gfns)
		sfns = cData{formatI(i)};
		sfns{1} = regexprep(sfns{1},'^\w*=','');
		if(i==1)
			mask = ~specialI;
		else
			specialI = find(specialI);
			mask = specialI(strcmp(specials,gfns{i}));
		end
		dat = vertcat(cData{mask});
		if(i>1)
			dat(:,1) = regexprep(dat(:,1),'\w*=','');
		end
		s.(gfns{i}) = cell2struct(dat,sfns,2);
	end
else
	if(sum(specialI)>sum(formatI)); warning('cymeread:unknownspecial','unknown special line in this section'); end
	% first line and some subsequent line are formats
	if(~isequal(formatI,specialI))
		warning('This table has multiple formats and a changing default format, which is something I haven''t been programmed to deal with yet');
		% What I have in mind is that you could modify the "with specials" case above to be smarter about how it deals with them so that you can specify a format for a given field more than once, and then pass in all of the cumulative header lines to the recursive call.  Or perhaps skip modifying the code above and instead be smarter about which header lines are passed in recursively.  But I don't even know how the format is supposed to work in such a case (e.g. which format becomes the default), so I'm not going to write the code yet.
	end
	headerlines = find(formatI);
	endlines = [headerlines(2:end)-1 length(formatI)];
	for i=1:length(headerlines)
		dat = vertcat(cData{headerlines(i)+1:endlines(i)});
		sfns = cData{headerlines(i)};
		sfns{1} = regexprep(sfns{1},'^\w*=','');
		sfns = sanitize(sfns);
		s{i} = cell2struct(dat,sfns,2);
	end
end

end

function fn = sanitize(fn)
fn = fnSanitize(fn,'x');
end
