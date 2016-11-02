%% Matlab Library
%
% Title: Histogram Bins
%
% Author: Bryan Urquhart
%
% Description:
%   This code will generate histogram bins based on a set of 2 bounds and a
%   number of bins. The left, right and central points of each bin are returned.

function [y i] = lmHistBinVal( bins , x )
%% Process Input Arguments

if( size(bins,2) ~= 3 )
  error( 'You must have left, center and right bins. See lmHistBins()' );
end

% Get number of bins
nbins = size(bins,1);

% Initialize the output
y = zeros( nbins , 1 );

%% Bin Data

i = nan(length(x),1);

% Loop over bins and determine how many elements lie within a given bin
for idx = 1:nbins
  % Count the number of data within a bin
  if( nargout == 2 )
    if( idx ~= 1 && idx ~= nbins )
      index = find( bins(idx,1) <= x & x < bins(idx,3) );
    else
      if( idx == 1 )
        index = find( x < bins(idx,3) );
      else
        index = find( bins(idx,1) <= x );
      end
    end
    i(index)=idx;
    y(idx) = numel(index);%sum( bins(idx,1) <= x & x < bins(idx,2) );
  else
    if( idx ~= 1 && idx ~= nbins )
      y(idx) = sum( bins(idx,1) <= x & x < bins(idx,3) );
    else
      if( idx == 1 )
        y(idx) = sum( x < bins(idx,3) );
      else
        y(idx) = sum( bins(idx,1) <= x );
      end
    end
  end
end
