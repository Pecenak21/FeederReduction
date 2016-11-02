function keymagic( shadowpath, powerimages, ptimes, starttime )
%Keymagic adds l/r/u/d arrow key functionality to the figure
timeindex = 1;
forecastindex = 1;
if ischar(starttime)
	if(length(starttime)==14 && all(starttime>='0') && all(starttime <= '9'))
		starttime = datevec(starttime,'yyyymmddHHMMSS');
	else
		starttime = datevec(starttime);
	end
elseif strcmp(class(starttime),'bu.util.Time')
	starttime = bu.util.Time.timeToDatevec(starttime);
	starttime = starttime';
elseif isnumeric(starttime)
	starttime = datevec(starttime);
else
	error('invalid start time');
end




h = figure;
myidx = ptimeindex(starttime);
dstr = sprintf('%04d%02d%02d',starttime(1),starttime(2),starttime(3));
% modyify these two to change which we're loading:
fdat = load(sprintf('%s/%s/cloudmap/cloudmap_%s%02d%02d%02d.mat',shadowpath,dstr,dstr,starttime(4),starttime(5),starttime(6)));
%fdat = fdat.cloudmap.rbrmap-fdat.cloudmap.cslmap;
fdat = siShadow(fdat.cloudmap,fdat.cloudmap.time,'map');
%sh = imagesc(fdat,'alphadata',(fdat~=1)); % you probably want this one for binary forecasts
sh = imagesc(fdat); % and this one for full-channel
hold on;
ph = imagesc(powerimages{myidx},'AlphaData',(powerimages{myidx}~=0)*0.7);

set(gca,'clim',[0 1]);

set(h,'keypressfcn',@doKeyDown);

	function doKeyDown(fig_h,eventdata)
		switch(eventdata.Key(1))
			case 'u'
				starttime = datevec(datenum(starttime + [0 0 0 0 10 00]));
			case 'd'
				starttime = datevec(datenum(starttime - [0 0 0 0 10 00]));
			case 'l'
				starttime = datevec(datenum(starttime - [0 0 0 0 0 30]));
			case 'r'
				starttime = datevec(datenum(starttime + [0 0 0 0 0 30]));
		end
		% show the new power data
		myidx = ptimeindex(starttime,myidx);
		set(ph,'cdata',powerimages{myidx},'alphadata',(powerimages{myidx}~=0)*0.7);
		% load the new cloud data
		dstr = sprintf('%04d%02d%02d',starttime(1),starttime(2),starttime(3));
		fdat = load(sprintf('%s/%s/cloudmap/cloudmap_%s%02d%02d%02d.mat',shadowpath,dstr,dstr,starttime(4),starttime(5),starttime(6)));
		fdat = siShadow(fdat.cloudmap,fdat.cloudmap.time,'map');
		%fdat = fdat.cloudmap.rbrmap-fdat.cloudmap.cslmap;
		% show the new cloud data
		set(sh,'cdata',fdat);
		fprintf('.');
	end

	function idx = ptimeindex(time,ihint)
		if(nargin < 2); ihint = 1; end;
		ttime = datevec(ptimes(ihint));
		while ~all(ttime == time)
			if datenum(time)>datenum(ttime)
				ihint = ihint+1;
			else
				ihint = ihint -1;
			end
			ttime = datevec(ptimes(ihint));
		end
		idx = ihint;
	end

end

