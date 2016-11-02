%% Matlab Visualization Library
%
%  Title: PNG to AVI
%
%  Author: Bryan Urquhart
%  Date  : September 26, 2012
%
%  Description:
%    Converts a folder filled with PNGs to an avi video file. A step can be
%    specified if not all pngs should be converted. If you only need specific
%    pngs, place them in a separate folder.
%  
function visPNGtoAVI( input , output , varargin )
%% Process Input Arguments

% Option defaults
o.step = 1;
o.overwrite = false;
o.p.fps = [];
o.p.compression = [];
o.p.quality = [];

% handle varargin
if( ~isempty( varargin ) )
  args = argHandler(varargin);
  for idx = 1:size(args,1)
    switch( args{idx,1} )
      case 'step'
        o.step = args{idx,2};
      case 'overwrite'
        o.overwrite = args{idx,2};
      case {'framerate','fps'}
        o.p.fps = args{idx,2};
      case 'compression'
        o.p.compression = args{idx,2};
      case 'quality'
        o.p.quality = args{idx,2};
    end
  end
end

% convert inputs
input  = java.io.File( char(input ) );
output = java.io.File( char(output) );

if( ~input.isDirectory() )
  error( 'Must specify valid input directory.' );
end

if( output.isFile() && ~o.overwrite )
  error( 'Output file exists and overwrite was not specified.' );
end

% Construct avifile param value array
params = {};
pnames = fieldnames( o.p );
for idx = 1:numel(pnames)
  if( ~isempty( o.p.(pnames{idx}) ) )
    params{end+1} = pnames{idx}; %#ok<*AGROW>
    params{end+1} = o.p.(pnames{idx});
  end
end


%% Load up directory

% List the files
file.list = input.listFiles();

% Find the pngs
tmp = java.util.Vector();
for idx = 1:file.list.length()
  if( ~file.list(idx).getName().endsWith('.png') ), continue; end
  tmp.add( file.list(idx) );
end
file.list = tmp.toArray( javaArray( 'java.io.File' , tmp.size() ) );
clear tmp;

%% Construct AVI

% set up the working avi object
attempt = false;
while(~attempt)
  try
    if( isempty(params) )
      aviobj = avifile(char(output)); %#ok<*TNMLP>
    else
      aviobj = avifile(char(output),params{:});
    end
    attempt = true;
  catch e %#ok<NASGU>
    % Clear open avi files
    clear mex;
    attempt = true;
  end
end
v = 1:o.step:file.list.length;
N = numel(v); clear v;
count = 0;

for idx = 1:o.step:file.list.length
  
  % User Message
  count = count + 1;
  fprintf( 1 , 'Adding frame %d of %d.\n' , count , N );
  
  % Extract the frame data from the image
  F.cdata = imread( char( file.list(idx)) );
  F.colormap = [];
  
  % Add to avi object
  aviobj = addframe(aviobj,F);
  
end

% Close the current avi object
aviobj = close(aviobj); %#ok<NASGU>

