%% Matlab Library
%  
%
%  Title: Image Scaling
%
%  Author: Bryan Urquhart
%
%  Description:
%    Generate a 3 band image from a single 2D matrix and a colormap. Returned
%    image has class uint8
%
function img = imagesc2imshow( img , cMap , varargin )
%% Process Input Arguments

% Verify that the image is 1 band (2d matrix)
if( length( size( img ) ) ~= 2 )
  error( 'Input image must be a single band' );
end

if( size( cMap , 2 ) ~= 3 )
  error( 'Colormap must be 3 columns corresponding to R,G & B' );
end

if( size( cMap , 1 ) < 2 )
  error( 'Colormap must contain at least 2 distinct colors' );
end

cLim = [];
if( ~isempty(varargin) )
  % Pass args to arg Handler
  args = argHandler(varargin);
  % Switch on args
  for i = 1:size(args,1)
    switch(lower(args{i,1}))
      case 'clim'
        cLim = args{i,2};
    end
  end
end

%% Scale the image

if( isempty( cLim ) )
  vMin = min(img(:));
  vMax = max(img(:));
else
  vMin = cLim(1);
  vMax = cLim(2);
end

N = size(cMap,1)-1;

img = uint8( (img - vMin)/(vMax - vMin) * N );
%img = ind2rgb( img , cMap );
img = uint8(ind2rgb( img , cMap ).*256);
