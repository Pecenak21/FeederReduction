%% Matlab Library
%
% Title: Sliding mean for linear time series
%
% Author: Bryan Urquhart
%
% Description:
%   This function generates a sliding mean for a linear data series. A
%   convolution with FFT is used for speed.
%
function y = lmFilter1( x , fSize , varargin )
%% Process Input Arguments

% Default filter size is central
type = 'central';
sigma = [];
shape = 'same';
x = x(:);
if( nargin == 1 )
	fSize = 10;
end

if( ~isempty( varargin ) )
  args = argHandler(varargin);
  for idx = 1:size(args,1)
    switch( args{idx,1} )
      case 'type'
        type = args{idx,2};
      case 'sigma'
        sigma = args{idx,2};
			case 'shape'
				shape = args{idx,2};
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
m = floor( numel(mfilt)/2 );
n = numel(x);

if(n<m), y = m; return; end

%% Apply filter

% Use filter2 to perform convolution
% Note we are only using the valid portion of the filtered signal so we replace
% all unfiltered values in the center of the signal with the originals values
% (this is thought to be a better practice than using zero padding).

switch( shape )
	
  case 'full'
		
	case 'same'
    y = filter2( mfilt , x );
    y(1:m) = x(1:m); y(n-m+1:n) = x(n-m+1:n);
    index = ~isnan(x) & isnan(y);
    y( index ) = x( index );
    
    
	case 'valid'
		n = numel(x);
		y = nan(n,1);
    y(m+1:n-m,1) = filter2( mfilt , x , 'valid' );
%     y(1:m) = nan;%x(1:m);
% 		y(n-m+1:n) = nan;%x(n-m+1:n);

  otherwise
    y = filter2( mfilt , x );
end
	


% % Put data back in that was removed in filtering
% index = ~isnan(x) & isnan(y);
% y( index ) = x( index );


