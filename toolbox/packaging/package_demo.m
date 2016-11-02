% script to prep a demo directory for shipment:
% * remove unneeded/proprietary files
% * pcode most of the others, except a few specific ones that we anticipate may need to be edited
% * zip up the whole package for distribution

% our new approach is to generate a list of all the files in the current directory that need to be added to a zip file.  That way no deleting is needed

%% check that git is clean
if(~isunix()), error('packaging:bad_os','This script only works on real computers\n'); end
[~, gitdiffs] = system('git status -s');
if(~isempty(gitdiffs))
	error('packaging:not_clean','Must be run on a clean git repository.  This includes no untracked files.');
end

%% read file lists
fprintf('Reading file lists\n');
removelist = {};
fid = fopen(which('remove.flist'));
while ~feof(fid)
	removelist{end+1} = fgetl(fid);
end
fclose(fid);

pcodelist = {};
fid = fopen(which('pcode.flist'));
while ~feof(fid)
	pcodelist{end+1} = fgetl(fid);
end
fclose(fid);

%% remove those files, along with our lists and ourselves
removelist = [removelist, {which('pcode.flist'), which('remove.flist'), mfilename('fullpath')}];

%% get the list of files
%  all files except for the ones on pcodelist and the ones we're deleting
fprintf('Gathering file lists\n');
% find sequence to skip paths we don't want:
skippatt = sprintf(' \\! -path "./%s*"',removelist{:});
% run find and split into cell arrays
[~, p] = system(['find . -iname "*.m"' skippatt]);
[~, o] = system(['find . \! -iname "*.m" \! -path "*\.git*" \! -type d' skippatt]);
p = regexprep(p, '^\./', '', 'lineanchors');
o = regexprep(o, '^\./', '', 'lineanchors');
p(end) = []; %remove the trailing newline before splitting
o(end) = [];
p = regexp(p,'\n','split');
o = regexp(o,'\n','split');
% move the files that we're not going to pcode from the 'to pcode' list to the 'other files' list
o = [o, p(ismember(p,pcodelist))];
p(ismember(p,pcodelist)) = [];
% add the .p files to the "other" list
o = [o, regexprep(p, '\.m$', '.p')];

%% do the pcode operation
%  We use 'pcode' to pcode all the files, and then go through and replace the mfiles containing code with mfiles containing just the documentation.  We also have to set the file date so we won't get warnings about the mfile being newer than the pfile.
fprintf('Generating pcodes\n');
pcode(p{:},'-inplace');
copyrightstring = ['% Copyright © Jan Kleissl''s Solar Resource Assessment and Forecasting Laboratory 2010-', datestr(now,'yyyy'), ' All Rights Reserved.'];
for i=1:length(p)
	fn = p{i};
	fdate = eval('java.io.File(fn).lastModified;');
	fhelp = help(fn);
	fhelp = regexprep(fhelp, '^(.)','%$1','lineanchors');
	fhelp = [fhelp sprintf('%%\n%s\n', copyrightstring)];
	fid = fopen(fn,'w');
	fprintf(fid, '%s', fhelp);
	fclose(fid);
	eval('java.io.File(fn).setLastModified(fdate);');
end

%% create the zip archive
fprintf('Creating zip archive\n');
% try to get a version number tag from git
[~, v] = system('git tag -l --contains HEAD');
if(isempty(v))
	v = 'vx_x';
else
	v = regexp(v,'\n','split'); v = v{1};
end
% create a filename and zip it up
fn = ['forecast_' v '.zip'];
zip(fn, [p, o]);

%% Cleaning up
fprintf('Cleaning up my mess!\n');
system('find . -iname "*.p" | xargs rm');
system('git checkout -- .');

fprintf('Done!\n');

