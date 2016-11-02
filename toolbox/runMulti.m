function runMulti(days, imager, deployment, rootdir, outputdir, notification_email, runstep)
% runMulti starts parallel historical forecasts for several days
%
% It does this by using repoMirror to create several mirrors of the code, and then tweaking several config parameters before starting the forecasts running in new sessions of 'screen'
%
% Inputs:
%	days - list of dates as array of MATLAB datenums or cellstring (required)
%	imager - the name of the imager to use (required)
%	deployment - the name of the deployment to use (optional, default is whatever's in global forecast.conf)
%	rootdir - location in which to create mirror repos and log files (optional, ~/run)
%	outputdir - different location to create output (default is rootdir)
%	notification_email - email to notify if the run fails.  if unspecified, no emails will be sent.  You may specify multiple emails by separating them with a comma.
%   runstep - step to run. default: 'all'
%
% If you need to specify a later input and want to use the default on an earlier one, just pass an empty value.
%
% Examples:
%		runMulti({'2012-11-11','2012-11-13'}, 'USI_1_2');
%		runMulti('2012-11-14', 'USI_1_2', '','', '/mnt/lab_18tb3/database/USI/analysis/newrun', 'solarucsd@gmail.com');
%
% Bugs:
%		currently we assume that forecasting is to be done in a time zone in which 08:00 UTC is during the night

%% Input processing
if(iscellstr(days) || ischar(days)), days = datenum(days); end
if(isa(imager,'siImager')), imager = imager.name; end
if(nargin < 3), deployment = ''; end
if(nargin < 4 || isempty(rootdir)), rootdir = '~/run'; end
if(nargin < 5 || isempty(outputdir)), outputdir = rootdir; end
if(nargin < 6), notification_email = ''; end
if(nargin < 7 || isempty(runstep)), runstep = 'all'; end
days = floor(days);
% check that we're in the repo root, assumed by the repoMirror script:
if(~exist('./def_addpath.m','file'))
	origd = cd(regexprep(which('def_addpath'),'def_addpath.m$',''));
end

%% Check status with git
%  We want to confirm if the user is not on the 'master' branch, and if the working copy is dirty, because these are states that are normally not used for large runs - usually the user will just run one day at a time

% What branch are we on?
[~, gitbranch] = system('git branch --list | grep ''\*''');
if(isempty(gitbranch) || ~isempty(regexp(gitbranch,'fatal: Not a git repository')))
	warning('runMulti:notGit', 'We are not currently working with code in a git repository.\nWorking with history-less code is not recommended because it makes results difficult to reproduce.');
	okaytorun = input('proceed even though code is not in git? y/N:', 's');
	if(isempty(regexpi(oktorun,'^\s*y'))), error('runMulti:notGit', 'canceled because code is not in git'); end
else
	if(isempty(gitbranch) || ~strcmp(gitbranch(3:end-1), 'master')) % if we're on the master branch, the above command should return '* master'
		warning('runMulti:notMaster', 'git reports that we are not on the master branch.  Current branch: %s', gitbranch(3:end-1));
		okaytorun = input('did you mean to forecast on a branch? Y/n:', 's');
		if(~isempty(regexp(oktorun,'^\s*n'))), error('runMulti:notMaster', 'canceled because we are not on the master branch'); end
	end

	% And whether the working copy is dirty
	%  Note that we ignore 'untracked' files here
	[~, gitdiffs] = system('git status -s -uno');
	if(~isempty(gitdiffs))
		warning('runMulti:notClean', 'git reports that the following files are changed.  Large runs without a clean working copy are not recommended!\n%s', gitdiffs);
		oktorun = input('proceed despite dirty working copy in git? y/N:','s');
		if(isempty(regexpi(oktorun,'^\s*y')))
			error('runMulti:notClean', 'canceled due to dirty working copy');
		end
	end
end

%% create a template string for sed to modify forecast.conf
% this will take care of imager, deployment, and email
outputdir = regexprep(outputdir, '/', '\\/'); % escape slashes
sedstr = ['-e ''s/^imager\s.*/imager		' imager '/'' -e ''s/^outputDir\s.*/outputDir		' outputdir '/'''];
if(~isempty(deployment))
	sedstr = [ sedstr ' -e ''s/^deployment\s.*/deployment		' deployment '/'''];
end
if(~isempty(runstep))
	sedstr = [ sedstr ' -e ''s/^step\s.*/step		' runstep '/'''];
end
if(~isempty(notification_email))
	sedstr = [ sedstr ' -e ''s/^emaillist\s[^#]*/emaillist		' notification_email '/'' -e ''s/sendemail\s[^#]*/sendemail		1/'''];
end

%% start each day
here = pwd;
for i=1:length(days)
	dayname = lower(datestr(days(i),'mmmdd'));
	logname = [rootdir '/' datestr(days(i),'yyyymmdd') '.log'];
	UTCoffset = -8/24;
	
	% create the new repo mirror and switch to it
	repoMirror([rootdir '/' dayname], true);
	cd([rootdir '/' dayname]);
	% modify forecast.conf for this day
	todaysed = [ '-e ''s/^startTime\s.*/startTime ' datestr(days(i)-UTCoffset, 31) '/'' -e ''s/^endTime\s.*/endTime ' datestr(days(i)-UTCoffset+1, 31) '/'''];
	system(['sed ' sedstr ' ' todaysed ' -i -- conf/local/forecast.conf']);
	% launch the forecast
	% screen flags:
	%	-d -m		detaches the new session so we aren't getting massively confused
	%	-S dayname	names the session properly
	% matlab flags:
	%	-nodisplay	non graphical
	%	-logfile xx	write session log
	%	-r xxx		run the command to start the forecast
	system(['screen -d -m -S ' dayname ' matlab -nodisplay -logfile ' logname ' -r ''def_addpath; siForecastMain''']);
	fprintf('started new matlab task to forecast for the day of %s\n', datestr(days(i), 31));
	% return to the main directory
	cd(here);
end

% return to the original directory if necessary
if(exist('origd','var')), cd(origd); end

end

