function repoMirror(newlocation, clean)
% repoMirror(newlocation, [clean]) makes a cheap mirror of the forecast code in a second location for running multiple configurations in parallel
%
%  This is accomplished by symlinking everything except for the local config files into the new location, (the local configs are copied) so that you can simply edit the local configs and run.  The symlink operation is essentially instantaneous
%  The optional 'clean' flag, can be used to tell repoMirror to use fresh copies of all config files (except piServer.conf) from the defaults instead of duplicating the local config.

% create the directory
if(exist(newlocation,'dir'))
	error('cannot overwrite existing location');
end
mkdir(newlocation);

% symlink everything
origd = cd(newlocation);
system(['ln -s ' origd '/* ./']);

% replace conf with a real directory
%delete conf % can't use the matlab delete because it resolves the symlink
system('rm conf');
mkdir conf
cd conf

% symlink all the default config files in
system(['ln -s ' origd '/conf/*.conf ./']);

% and copy all the files for the local conf
mkdir local
cd local
if(exist('clean','var') && clean)
	% copy default conf
	system(['cp -H ' origd '/conf/*.conf ./']);
	% and local pi server conf
	system(['cp -H ' origd '/conf/local/piServer.conf ./']);
else
	% just copy the local conf
	system(['cp ' origd '/conf/local/*.conf ./'])
end

% finally, cd back to the original location:
cd(origd);

end
