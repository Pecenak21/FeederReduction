function r = rbr(img, offset)
% RBR is the Red-Blue-Ratio of an image, optionally with a zero offset
% subtracted from data before taking the ratio. May 24, 2013: Dark offset
% now subtracted in siImager. No offset will be subtracted regardless of
% second input to prevent duplicate offset subtraction.

offset = 0;

if(nargin > 1) % offset specified
	r = (double(img(:,:,1))-offset)./(double(img(:,:,3))-offset);
	r(img(:,:,3)==0) = nan;
else
	r = double(img(:,:,1))./double(img(:,:,3));
end