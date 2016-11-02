%% Matlab Library
%
% Title: Sliding std dev for linear time series
%
% Author: Bryan Urquhart
%
% Description:
%   This function generates a sliding std dev for a linear data series. A
%   convolution with FFT is used for speed.
%
function y = lmDev1( x , fSize , varargin )
%% Process Input Arguments

% Default filter size is central
type = 'central';
sigma = [];

if( ~isempty( varargin ) )
  args = argHandler(varargin);
  for idx = 1:length(args)
    switch( args{idx,1} )
      case 'type'
        type = args{idx,2};
      case 'sigma'
        sigma = args{idx,2};
    end
  end
end

%% Create filter

% Create kernel
switch( type )
  case 'central'
    nfilt = floor( fSize / 2 );
    mfilt = [ ones(nfilt,1) ; 1 ;  ones(nfilt,1) ];
  case 'lagging'
    mfilt = [ ones(fSize-1,1) ; 1 ; zeros(fSize-1,1) ];
  case 'leading'
    mfilt = [ zeros(fSize-1,1) ; 1 ; ones(fSize-1,1) ];
  case 'gaussian'
    
    if( isempty(sigma) )
      sigma = 0.5;
    else
      if( sigma < 0 )
        error( 'lmFilter: sigma must be greater than zero.' );
      end
    end
    
    if( mod(fSize,2) )
      X = (ceil(-fSize/2):floor(fSize/2))';
    else
      X = (-(fSize-1)/2:1:(fSize-1)/2)';
    end
    
    a = 1/( sigma * 2 * pi );
    b = - X.^2 / ( 2 * sigma^2 );
    mfilt = a*exp(b);
    
  otherwise
    nfilt = floor( fSize / 2 );
    mfilt = [ ones(nfilt,1) ; 1 ;  ones(nfilt,1) ];
end

% Ensure that we don't change signal power content - normalize kernel to unity
mfilt = mfilt / sum(mfilt(:));

%% Apply filter

% Compute the local mean
x_ = lmFilter1( x , fSize , varargin{:} );

% Compute the local mean square deviation
y  = lmFilter1( (x - x_).^2 , fSize , varargin{:} );

% Take the square root
y = sqrt( y );


