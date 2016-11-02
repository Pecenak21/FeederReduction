%% Matlab Library
%  Bryan Urquhart
% 
%  Convert avi into gif animation
function avi2gif( input , output )

% Get multimedia object
obj = mmreader(input);

% Verify object is valid
if( ~obj.isvalid )
  error( 'Input file does not yield valid mmreader object!' );
end

% Get file type
filetype = obj.Name(end-2:end);
switch( filetype )
  case 'wmv'
    img.nFrames = round( obj.Duration * obj.FrameRate );
    img.nFramesValid = false;
  otherwise
    img.nFrames = obj.NumberOfFrames;
    img.nFramesValid = true;
end

% Get info from object
img.height  = obj.Height;
img.width   = obj.Width;

% Set the time for each frame to consume
gif.delaytime = .125/obj.FrameRate; % half a second per frame

% Read one frame at a time.
h_wait = waitbar(0,'Converting AVI to GIF - Please wait...');
idx = 1;

while( true ) %for idx = 1:img.nFrames

  % Get a frame from the movie 
  % Handle exceptions that occur from reading too many frames (wmv files)
  try %#ok<ALIGN> 
    img.cdata = obj.read( idx );
  catch e, break; end %#ok<NASGU>
  
%   close 1;
%   figure(1); imshow(img.cdata); pause(0.1); disp(idx);
%   idx = idx + 5;
  
  % Create indexed gif image and map
  [gif.img,gif.map] = rgb2ind(img.cdata,256);
 
  % Write a gif frame to the output file
  if( idx == 1 )
    imwrite(gif.img,gif.map,char(output),'gif', ...
                                          'DelayTime' , gif.delaytime , ...
                                          'WriteMode' , 'overwrite' );
  else
    imwrite(gif.img,gif.map,char(output),'gif', ...
                                          'DelayTime' , gif.delaytime , ...
                                          'WriteMode' , 'append' );
  end
  
  % Update waitbar
  waitbar( idx / img.nFrames , h_wait , ['Converting AVI to GIF - Frame ' num2str(idx)]);
  
  % Increment index
  idx = idx + 2;
  
  % Check for exit condition
  if( img.nFramesValid )
    if( idx > img.nFrames ), break; end
  end
end
% Close the waitbar
close(h_wait);