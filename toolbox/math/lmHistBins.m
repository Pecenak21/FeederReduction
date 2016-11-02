%% Matlab Library
%
% Title: Histogram Bins
%
% Author: Bryan Urquhart
%
% Description:
%   This code will generate histogram bins based on a set of 2 bounds and a
%   number of bins. The left, right and central points of each bin are returned.

function X = lmHistBins( nbins , bounds )
%% Process Input Arguments

if( rem( nbins , 1 ) )
  error( 'The number of bins must be an integer' );
end

if( length(bounds) ~= 2 )
  error( 'You can only specify two boundary points' );
end

if( bounds(1) > bounds(2) )
  error( 'Bounds must be increasing' );
end

% Initialize the output
X = [];

%% Generate bins

range = bounds(2) - bounds(1);
binsize = range / nbins;

left   = (bounds(1):binsize:bounds(2)-binsize)';
center = left + binsize/2;
right  = left+binsize;

X = [left center right];
