%% Matlab Library
%
%  Title: Image Overlay
%
%  Author: Bryan Urquhart
%
%  Description:
%    To evaluate the maps generated for sky images there was a need to overlay a
%    data matrix on top of a sky image. This function generalizes the process of
%    overlaying data from matrix onto another (which can be an image).
%
function hf = imageoverlay( img , overlay , varargin )
%% Process Input Arguments

% Make sure that overlay is only 2D
if( numel( size( overlay ) ) ~= 2 )
  error( 'Overlay can only be two dimensional' );
end


% Validate the sizes of the input images
if( numel(size(img)) > 3 )
  error( 'You can not plot an image with more than 3 dimensions' )
end
if( numel(size(img)) == 3 )
  [M N ~] = size(img);
else
  [M N] = size(img);
end
if( size( overlay , 1 ) ~= M || size( overlay , 2 ) ~= N )
  error( 'The number of rows and columns in the image and the overlay must be identical.' );
end


CLim = [ min(overlay(:)) max(overlay(:)) ];
opacity = 0.5;
CMapImg = gray(256);
CMapOverlay = jet(256);

if( ~isempty(varargin) )
  % Pass args to arg Handler
  args = argHandler(varargin);
  % Switch on args
  for i = 1:size(args,1)
    switch(lower(args{i,1}))
      case 'clim'
        CLim = args{i,2};
      case 'opacity'
        opacity = args{i,2};
      case 'cmapimg'
        CMapImg = args{i,2};
      case 'cmapoverlay'
        CMapOverlay = args{i,2};
    end
  end
end

%% Store opengl status
glstatus = opengl('data');
if( ~glstatus.Software )
  opengl('software');
  disp('Setting opengl to software mode' );
end

%% Visualization

% Generate transparency map
imAlphaData = repmat(opacity,[M N]);
imAlphaData(isnan(overlay)) = 0;


% Some example code
hf = figure(gcf);
set( hf , 'units' , 'normalized' );
set( hf , 'position' , [.2 .2 .6 .6] );
set( hf , 'color' , [ 1 1 1 ] );

% Now set up axes that overlay the background with the image
% Notice how the image is resized from specifying the spatial 
% coordinates to locate it in the axes.

imagesc(overlay);
axis image;
hold on;
imagesc(img); colormap( CMapImg );
% Overlay the image, and set the transparency previously calculated
h2 = imagesc(overlay); colormap(CMapOverlay); colorbar; caxis(CLim); 
set(h2,'AlphaData',imAlphaData);


%% Reset OpenGL status
pause(0.001);
if( ~glstatus.Software )
  opengl('hardware');
  disp('Setting opengl to hardware mode' );
end
