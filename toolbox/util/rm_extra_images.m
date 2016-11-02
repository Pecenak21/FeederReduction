function rm_extra_images(imager, day)
% remove images stored before sunrise or after sunset for a given imager
%
% images are actually moved to a folder called pending_delete in the parent folder of the images for that date

error('you cannot use this function unless you really really want to!');
% If you want to remove extra images, you need to comment out the line
% above.  This is to prevent accidental loss of data.
%
% This function comes with no guarantees that it will not harm you or your
% children or your advisor's children

% get day in matlab datenum format for the start of the day (in UTC)
if(ischar(day) || iscell(day))
	day = datenum(day);
end
day = floor(day);


% calculate the solar angles all day each day, and then compare to get 
for i = 1:length(day);
	day_ = day(i);
	imgdir = [imager.imageDir(day_) '/' datestr(day_,'yyyymmdd')];
	if(~exist([imager.imageDir(day_) '/pending_delete'],'dir'))
		mkdir([imager.imageDir(day_) '/pending_delete']);
	end
	owd = cd(imgdir);
	flist = dir('2*.*');
	ft = datenum({flist.name},'yyyymmddHHMMSS');
	p = siSunPosition(ft, imager.position);
	isNight = [p.zenith]>93;
	
	fx = fopen('/tmp/files2delete.txt','w');
	fprintf(fx,'%s\n',flist(isNight).name);
	fclose(fx);
	
	system('cat /tmp/files2delete.txt | xargs -I {} mv {} ../pending_delete')

	cd(owd);
end
