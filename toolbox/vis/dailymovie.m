function dailymovie(imager,day)
% DAILYMOVIE(imager,day)
% generate a daily movie for the specified imager on the given day
%	makes use of the dailymovie.pl script and x264.  This more or less makes this a non-windows function
%
% input:
%			day : either datenum format or '2012-11-11 11:11:11' format.

% we need a day as a datenum because we need to give the perl script the correct starting hour
if(nargin < 2)
	% default to today in whatever the local timezone is; this is a bit sketchy
	day = now;
elseif(ischar(day))
	day = datenum(day);
end
% start at the beginning of a day
day = floor(day);

% name the output file
output = getfield(readConf(siGetConfPath('setup.conf'),0),'MOVIE_DIR'); %#ok<GFLD>
output = siNormalizePath(output);
iname = regexprep(imager.name,'[Uu][Ss][Ii]_?(\d+)_(\d+)','usi$1-$2');
output = sprintf('%s/%s',output,iname);
if(~exist(output,'dir'))
	mkdir(output);
end
output = sprintf('%s/%s.mp4',output,datestr(day,'yyyymmdd'));
if(exist(output,'file'))
	error('The automatically chosen output file name already exists.  Cancelling');
end

% adjust for UTC time on the day (15 degrees of longitude per hour; just need to get close enough that this is in the middle of the night)
day = day - (imager.longitude/15)/24;

% run the perl script
fprintf('Preparing frames...\n');
[framedir] = perl('dailymovie_helper.pl', imager.imageDir(day), datestr(day,31), output);

% encode the video
if(~hasx264())
	error('dailymovie:missingX264','dailymovie requires x264 to be installed.  try ''sudo apt-get install x264''.\nYou will also need to do something like ''ln -s /usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.16 /usr/local/MATLAB/R2012b/sys/os/glnxa64/libstdc++.so.6''');
end
unix(['x264 ' framedir '/%04d.jpg --crf 18 -o ' output]);


end

function f = hasx264()
% checks for presence of x264
[~, p] = unix('which x264');
f = ~isempty(p);

end
