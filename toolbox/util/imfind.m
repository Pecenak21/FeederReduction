function [offset, max_cc, center] = imfind(img, pattern)
% IMFIND searches an image for the region that best matches a given pattern
%
% It acts as a simple wrapper around normxcorr2, similar to siNCC.
% 
% IMFIND returns an offset - the indices of the top left corner of the matched region
%	the max cross correlation - closer to 1 the better the match is
%	and the coordinates of the center of the matched region (not guaranteed to be an integer)

%% compute normalized cross correlation and find the max
cc = normxcorr2(pattern, img);
[max_cc, imax] = max(abs(cc(:)));
[row, col] = ind2sub(size(cc),imax(1));

%% correct to get the location as specified.
offset = [row, col] - size(pattern)+1;
center = offset + (size(pattern)-1)/2;

end
