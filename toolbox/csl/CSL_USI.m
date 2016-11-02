function [CSL path] = CSL_USI(imager, startTime, endTime, path)
% CSL_USI(imager, startTime, endTime, path) generates a new CSL for the imager using the specified clear day
%	All parameters are optional and will be inferred from forecast.conf if not specified.
%
%	The generated CSL will have automatically specified cutoffs that are
%	chosen in such a way as to match the previous CSL for that day (or for
%	a prior day) if possible
%
%	imager, startTime, and endTime can be specified as strings, e.g.
%	CSL_USI('USI_1_2', '2013-01-30 08:00:00', '2013-01-31 08:00:00')
%
%	you can also specify a batch of start times and end times to use clear images from multiple time periods:
%	CSL_USI('USI_1_2', {'2013-01-30 08:00', '2013-05-13 08:00'}, {'2013-01-31 08:00', '2013-05-14 08:00'})
%
%	If only a start time is specified, endTime will be inferred to match
%	startTime (thus using _only_ the images taken at exactly the times
%	specified as startTime), rather than read from forecast.conf, as would
%	be true if no startTime is specified either.
%

%% Initialize Variables
conf = readConf(siGetConfPath('forecast.conf'));

CSLRES = 1;

% Load the requested imager(s)
if(nargin < 1 || isempty(imager))
	imager = siImager(conf.imager);
elseif(ischar(imager))
	imager = siImager(imager);
end

% Times are UTC
if(nargin < 3 || isempty(endTime))
	if(nargin < 2 || isempty(startTime))
		% no times specified, infer everything from forecast.conf
		startTime = conf.startTime;
		endTime = conf.endTime;
	else
		% start time specified, assume we want a csl built from only the specified times
		endTime = startTime;
	end
end
if(ischar(startTime) || iscellstr(startTime)); startTime = datenum(startTime); end
if(ischar(endTime) || iscellstr(endTime)); endTime = datenum(endTime); end
time.begin = startTime;
time.end = endTime;
time.now = [];

% output path
if(nargin < 4 || isempty(path))
	dpath = [imager.confDir '/csl'];
	if( ~exist(dpath,'dir') )
		mkdir(dpath);
	end
	path = sprintf('%s/%s_CSL_a.mat', dpath, datestr(time.begin(1), 'yyyymmdd'));
	pathchar = 'a';
	while(exist(path,'file'))
		pathchar = char(pathchar+1);
		if(pathchar=='z'), error('out of automatic filenames to save under'); end
		path = sprintf('%s/%s_CSL_%c.mat', dpath, datestr(time.begin(1), 'yyyymmdd'),pathchar);
	end
end


%% Generate CSL
% allocate a struct to hold temp tables while we work
csl_t = struct('rbr',cell(conf.solarZenithMax,2),'rgb',[],'num',[]);
prt = nan;
for idx_time=1:numel(time.begin)
	fprintf('running from %s to %s\n',datestr(time.begin(idx_time)),datestr(time.end(idx_time)));
	
	% get a first image for this time
	% note that if we already found one from a previous time entry and it's equal to the current begin time, we just use that instead
	if(isempty(time.now) || time.now~=time.begin(idx_time))
		[~, time.now] = imager.nextImage(time.begin(idx_time)-1.5/24/3600, max(endTime)); % start searching images 1.5sec before specified begin time, so that we definitely get the image at that time, even if it might have been captured a second early
	end
	% update geometry mapping and zenith angles.  Assume these don't change within one time block
	image_zenith = imager.geomMap(time.now);
	mask = image_zenith <= 90; %conf.solarZenithMax;
	
	while(time.now <= time.end(idx_time))
		% Check if the sun is in a part of the sky where we can reasonably forecast
		sun = siSunPosition(time.now, imager.latitude, imager.longitude, imager.altitude);
		
		if( sun.zenith > conf.solarZenithMax )
			fprintf('Skip due to small solar elevation angle (%.2f) at time: %s\n', 90-sun.zenith , datestr(time.now));
			[~, time.now] = imager.nextImage(time.now);
			continue;
		end
		
		% Now that we're sure we need to do things, read the image from disk and crop it
		fprintf('Process at: %s. Zenith angle: %.f\t', datestr(time.now, 'yyyy-mm-dd HH:MM:SS'), sun.zenith);
		timer_ = tic;
		img1 = imager.readImage(time.now);
		
		% separate out the channels of the image
		RBR = rbr(img1,imager.darkoffset);
		
		% sun-pixel angle needed for coordinate remapping and blooming stripe search
		alpha = imager.sunPixelAngles( sun.zenith, sun.azimuth );
		
		% find and remove the blooming stripe
		% we only remove it in the RBR for now because we use the RBR to generate the mask below, so doing this should copy over to the RGB images as well.
		[mm] = siFindBloomingStripe(img1, alpha);
		if(~isempty(mm))
			RBR(:,mm) = nan;
		end
		% proceed with mapping into image coordinates
		mm = isnan(RBR) | isinf(RBR);
		ii = find(~mm & mask);
		I = sub2ind([90*CSLRES+1 180*CSLRES+1], round(image_zenith(ii)*CSLRES)+1, round(alpha(ii)*CSLRES)+1);
		RBR(mm) = 0;
		
		% Look up the temp table to use for the current image
		SZA = round(sun.zenith);
		prt = 1+(sun.azimuth > 180);
		if( isempty(csl_t(SZA,prt).rbr) )
			% create a new table
			RGB_temp = zeros((90*CSLRES+1)*(180*CSLRES+1), 3);
			RBR_temp = zeros(90*CSLRES+1, 180*CSLRES+1);
			num_rbr = RBR_temp;
		else
			% pull the temp table back out; we do this once to avoid extra indexing in the loop
			RGB_temp = csl_t(SZA,prt).rgb;
			RBR_temp = csl_t(SZA,prt).rbr;
			num_rbr = csl_t(SZA,prt).num;
		end
		
		% save the image data for the current step
		[y, z] = grouping_unique(I);
		img1 = reshape(img1,[size(img1,1)*size(img1,2) 3]);
		zn = 1;
		for i = y(:)';
			x = ii(z{zn});
			RBR_temp(i) = RBR_temp(i) + sum(RBR(x));
			num_rbr(i) = num_rbr(i) + numel(x);
			RGB_temp(i,:) = RGB_temp(i,:) + sum(double(img1(x,:)));
			zn = zn+1;
		end
		
		% Store the temp data back in the in-progress CSL
		csl_t(SZA,prt).rgb = RGB_temp;
		csl_t(SZA,prt).rbr = RBR_temp;
		csl_t(SZA,prt).num = num_rbr;
		
		% User message and next loop
		fprintf('Time elapsed: %.2f seconds.\n', toc(timer_));
		[~, time.now] = imager.nextImage(time.now, max(endTime));
	end
end
% quick check for no images is that prt is nan
if(isnan(prt))
	error('siCSL:noImages','There don''t appear to be any images in your selected time range');
end

%% Average the temp tables and save the result into the final structure
for prt=1:2
	for SZA=1:conf.solarZenithMax
		num_rbr = csl_t(SZA,prt).num;
		if(isempty(num_rbr)), continue; end
		s = size(num_rbr);
		CSL.prt(prt).SZA(SZA).Red = reshape(csl_t(SZA,prt).rgb(:,1),s)./num_rbr;
		CSL.prt(prt).SZA(SZA).Green = reshape(csl_t(SZA,prt).rgb(:,2),s)./num_rbr;
		CSL.prt(prt).SZA(SZA).Blue = reshape(csl_t(SZA,prt).rgb(:,3),s)./num_rbr;
		CSL.prt(prt).SZA(SZA).RBR = csl_t(SZA,prt).rbr./num_rbr;
	end
end

%% Set parameters and save
% some default cutoffs
CSL.cutoffs = [0.3 0.301];
% try to get a better set from a recent CSL
if(~isempty(imager.csl))
	CSL.cutoffs = imager.csl(find(([imager.csl.time] - startTime(1))<0, 1, 'last')).cutoffs;
end
% imager name and today's date
CSL.imager = imager.name;
CSL.generated = datestr(nowUTC, 31);
CSL.starttime = datestr(startTime, 31);
CSL.endtime = datestr(endTime, 31);

save(path,'-struct','CSL')

end
