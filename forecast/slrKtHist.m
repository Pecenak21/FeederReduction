%% Solar Resource Assessment Library
%
%  Title: Plot ramp rates
%
%  Author: Bryan Urquhart
%
%  Description:
%    This code will generate a histogram of ramp rates
function ktHist = slrKtHist( kt , varargin )
%% Process Input Arguments

% Set the default number of bins
nbins = 100;

% % Set the minimum threshold. 1 is nominally equal to 1 W or 1 W/m2
% minThresh = 0;


if( ~isempty( varargin ) )
  args = argHandler(varargin);
  for idx = 1:size(args,1)
    switch( args{idx,1} )
%       case 'threshold'
%         minThresh = args{idx,2};
      case 'numberofbins'
        nbins = args{idx,2};
    end
  end
end

%% Bin data


% Get the reduced dataset
% val = ramprate{idx}( abs(ramprate{idx}) > minThresh ) / pnom(idx);

% Get the min and max ramps
bmin = min( kt(:) );
bmax = max( kt(:) );

% Set up the bins
bins = lmHistBins( nbins , [bmin bmax] );

% Determine the number of elements in each bin
binval = lmHistBinVal( bins , kt(:) );

ktHist.bins   = bins;
ktHist.binval = binval;



