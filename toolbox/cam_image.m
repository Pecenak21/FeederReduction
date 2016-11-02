classdef cam_image < handle
	%CAM_IMAGE is a helper class to do things with images from the USI
	
	properties
		datetime;
		etime;
		imgdata;
		note = '';
		center = [];
		radius = NaN;
		dark = 38.15;
		bitdepth = 12;
	end
	
	properties (SetAccess = protected)
		hist_cache = [];
		hist_cache_bych = [];
		hist_is_crop = [false false false false];
	end
	
	methods
		function im = cam_image(time_str,etime)
			im.etime = etime;
			im.datetime = time_str;
			im.imgdata = imread(sprintf('~/usi_data/%s/%s_%04d.png',time_str(1:8),time_str,etime));
			if (strcmp(class(im.imgdata),'uint16'))
				im.imgdata = im.imgdata/16;
			else
				im.dark = im.dark/16;
				im.bitdepth = 8;
			end
			im.imgdata = double(im.imgdata);
		end
		
		function ri = rbr(im)
			if( im.is_cropped() )
				idat = im.cropped();
			else
				idat = im.imgdata;
			end
			ri = (idat(:,:,1)-im.dark)./(idat(:,:,3)-im.dark);
		end
		
		function fl = is_cropped(im)
			fl = ~(isempty(im.center) || length(im.center)~=2 || isnan(im.radius));
		end
		function ci = cropped(im)
			if( ~im.is_cropped() )
				warning('cam_image:nosize','Tried to crop sky image with missing center or radius');
				ci = im.imgdata;
			else
				ci = sky_image_crop(im.imgdata,struct('x',im.center(1),'y',im.center(2)),im.radius);
			end
		end
		function hi = hist(im,chnum)
			is_crop = im.is_cropped();
			if(nargin<2)
				chnum=0;
			end
			% check if we can use a cached histogram
			if(im.hist_is_crop(chnum+1) == is_crop)
				if chnum==0 && ~isempty(im.hist_cache)
					hi = im.hist_cache;
					return;
				elseif ~isempty(im.hist_cache_bych) && ~any(isnan(im.hist_cache_bych(:,chnum)))
					hi = im.hist_cache_bych(:,chnum);
					return;
				end
			end
			im.hist_is_crop(chnum+1) = is_crop;
			
			% get the image data
			if( is_crop )
				idat = im.cropped();
			else
				idat = im.imgdata;
			end
			if(chnum ~= 0)
				idat = idat(:,:,chnum);
			end
			
			% calculate the histogram
			hi = hist(idat(:),0:4095);
			
			% cache it for next time
			if chnum == 0
				im.hist_cache = hi;
			else
				if isempty(im.hist_cache_bych)
					im.hist_cache_bych = nan(length(hi),3);
				end
				im.hist_cache_bych(:,chnum) = hi;
			end
		end
		function rh = rbr_hist(im)
			d = im.rbr();
			[rh,b] = hist(d(:),2000);
			rh = [rh(:),b(:)];
		end
		
		function sp = sunstripe(im)
			sp = find(im.imgdata(1,:,3)>500);
		end
		function flag = has_sunstripe(im)
			flag = ~isempty(im.sunstripe());
		end
		
		function autocrop(im)
			% start by finding the edges in the x direction
			m = mean(im.imgdata(1000:1048,:,3));
			x = max(m); if x<180; x = 180; end;
			a = find(m(1:500)<x*0.25,1,'last');
			b = find(m(600:end)<x*.25,1)+600;
			im.radius = (b-a)/2;
			im.center(1) = (a+b)/2;
			% then the edges in the y direction
			m = mean(im.imgdata(:,1000:1048,3),2);
			x = max(m);
			a = find(m(1:500)<x*0.25,1,'last');
			b = find(m(600:end)<x*.25,1)+600;
			im.center(2) = (a+b)/2;
			disp(im.radius);
			disp((b-a)/2);
			im.radius = round(mean([im.radius, (b-a)/2]));
			im.center = round(im.center);
		end
	end
	
end

