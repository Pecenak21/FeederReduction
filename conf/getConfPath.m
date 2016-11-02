function lConfFile = getConfPath(cfName)
% getConfPath looks up the full path name to the configuration file of the requested name
% More importantly, getConfPath has various features to help the user maintain the relationship between the default config files and their local config files.  In particular:
%	* if there is no local conf file, the default will be copied there, and the user will be prompted whether they would like to edit the file
%	* if the default file is newer than the local file, the user is prompted if they would like to view a diff of the two files.  This is done using the external utility 'meld' on linux systems on which it is installed (since meld allows the user to edit the files), and with matlab's built-in visdiff otherwise.  Note that visdiff will fail in -nodisplay mode.
%
% siGetConfPath allows you to specify the conf file name with or without the .conf extension.  So both
%	getConfPath('sim')
%	getConfPath('sim.conf')
% return
%	'path_to_gridIntegration/conf/local/sim.conf'

% get a path to the conf directory
basep = normalizePath('$$');
basep = [basep '/conf'];
if( ~exist([basep '/local'],'dir') )
	mkdir([basep '/local']);
end

% determine if we can locate the conf file
dConfFile = [basep '/' cfName];
lConfFile = [basep '/local/' cfName];
% if we really can't find the file, try appending .conf to the end to see if that helps
if( ~exist(dConfFile,'file') && ~exist(lConfFile,'file') )
	dConfFile = [dConfFile '.conf'];
	lConfFile = [lConfFile '.conf'];
	if( ~exist(dConfFile,'file') && ~exist(lConfFile,'file') )
		fprintf('Looked for conf files named %s and %s\n',dConfFile,lConfFile);
		error('siGetConfPath:badConfName','No conf file of the requested name ''%s'' exists.',cfName);
	end
end

% we'll use this for some of our warnings. Warnings with html flags in them
% are a nice way to warn GUI users and give them several easily selectable
% options, but they don't work so well for CLI users, so we'll give _them_
% an input prompt about what to do, or at least a warning without all the
% HTML gunk in it.
isGUI = feature('ShowFigureWindows');

% now we have the paths we want, so proceed with checking that the local one exists
if ~exist(lConfFile,'file')
	copyfile(dConfFile,lConfFile);
	warning('siGetConfPath:usingDefaults','The local config file ''%s'' was missing and has been created using the default.  There is a good chance that the default will not work as desired.',cfName);
	opt = input('\n Would you like to:\n 1) Proceed anyway (using the defaults)\n 2) Stop and edit the file\n 3) just Stop.\n[1] ');
	if( opt == 3 )
		error('siGetConfPath:stop','User canceled');
	elseif( opt == 2)
		edit(lConfFile);
		error('siGetConfPath:stop_edit','User canceled to edit conf file');
  end
  lConfFile = java.io.File(lConfFile);
  lConfFile = char(lConfFile.getPath());
	return;
end

% finally check that the local file is at least as new as the default.
if ~exist(dConfFile,'file') % in case we got here with just a local and no default; should be rare
	warning('siGetConfPath:noDefault','You are using a conf file for which there is no default.  Please create a default file to help out other users!');
	return;
end
dT = dir(dConfFile); dT = dT.datenum;
lT = dir(lConfFile); lT = lT.datenum;
if(lT < dT)
	warningid = 'siGetConfPath:updatedDefaults';
	mesg = 'The local version of the ''%s'' conf file is older than the default. Some functions may not work if new settings have been added.';
	if(isGUI)
		warning(warningid,[mesg ' (<a href="matlab:visdiff(''%s'',''%s'');">view a diff</a> | <a href="matlab:java.io.File(''%s'').setLastModified(java.lang.System.currentTimeMillis);">update the file date</a> | <a href="matlab:warning(''off'',''%s'');">disable this warning</a>)\n'], cfName, dConfFile, lConfFile, lConfFile, warningid);
		opt = warning('query',warningid);
		if(strcmp(opt.state,'on') )
			pause(2);
		end
	else
		warning(warningid,[mesg '  To suppress this message, you can run "touch %s" in a system shell\n'],cfName,lConfFile);
	end
end

lConfFile = java.io.File(lConfFile);
lConfFile = char(lConfFile.getPath());

end