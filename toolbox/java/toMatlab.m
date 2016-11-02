%% Matlab Library
%  Java subset
%  Author: Bryan Urquhart
%  Date:   July 12, 2011
%
%  Description:
%    Converts a java style array to a matlab style array
%
function mArray = toMatlab( jArray , varargin )
%% Input process and validation

% Initialize the return param to an empty array
mArray = [];

% Just exit if the array is empty
if( isempty( jArray ) ), return; end;

% Verify the incoming data structure is a java structure
if( ~isjava(jArray) )
  error( 'Input argument is not a java array!' );
end

% Get the array sizing
N = size (jArray);
P = numel(jArray);

% Make javaArray into single column
jArray = jArray(:);

% If the array hasn't been populated the colon operator will return an empty
% array. We will just throw an error in this case
if( isempty( jArray ) )
  error('The java array input was empty!');
end

% Check if the user specified a particular conversion
type = '';  
typeCheck = false;
if( ~isempty( varargin ) )
  
  % Extract the first argument
  arg = varargin{1};
  
  % Convert to a char if it's a java string
  if( isa(arg,'java.lang.String') ), arg = char(arg); end
  
  % Check to make sure the argument is a char before proceeding
  if( ~isa(arg,'char') )
    error( 'Input argument must be a ''char''.' );
  end
    
  % Set the type to the input argument
  type = arg;
  
  % Set type checking flag
  typeCheck = true;
end

%% Perform cast to matlab type

% Allocate matlab array
mArray = nan( N );

% Cast to desired type. If the user doesn't specify a type, just assume double
if( ~typeCheck )
  
  try
    for idx = 1:P
      mArray(idx) = jArray(idx).doubleValue();
    end
    mArray = double( mArray );
  catch e
    mArray = nan( N );
    for idx = 1:P
      mArray(idx) = double(jArray(idx));
    end
    mArray = double( mArray );
  end
else
  
  % Switch on type
  switch( type )

    % Logical array of true and false values
    case 'logical'
      for idx = 1:P
        mArray(idx) = jArray(idx).doubleValue();
      end
      mArray = logical(mArray);

    % Characters array
    case 'char'
      for idx = 1:P
        mArray(idx) = jArray(idx).shortValue();
      end
      mArray = char(mArray);

    % 8-bit signed integer array
    case 'int8'
      for idx = 1:P
        mArray(idx) = jArray(idx).intValue();
      end
      mArray = int8(mArray);

    % 8-bit unsigned integer array
    case 'uint8'
      for idx = 1:P
        mArray(idx) = jArray(idx).intValue();
      end
      mArray = uint8(mArray);

    % 16-bit signed integer array
    case 'int16'
      for idx = 1:P
        mArray(idx) = jArray(idx).intValue();
      end
      mArray = int16(mArray);

    % 16-bit unsigned integer array
    case 'uint16'
      for idx = 1:P
        mArray(idx) = jArray(idx).intValue();
      end
      mArray = uint16(mArray);

    % 32-bit signed integer array
    case 'int32'
      for idx = 1:P
        mArray(idx) = jArray(idx).intValue();
      end
      mArray = int32(mArray);

    % 32-bit unsigned integer array
    case 'uint32'
      for idx = 1:P
        mArray(idx) = jArray(idx).intValue();
      end
      mArray = uint32(mArray);

    % 64-bit signed integer array
    case 'int64'
      for idx = 1:P
        mArray(idx) = jArray(idx).longValue();
      end
      mArray = int64(mArray);

    % 64-bit unsigned integer array
    case 'uint64'
      for idx = 1:P
        mArray(idx) = jArray(idx).longValue();
      end
      mArray = uint64(mArray);

    % Single-precision floating-point array
    case 'single'
      for idx = 1:P
        mArray(idx) = jArray(idx).floatValue();
      end
      mArray = single(mArray);

    % Double-precision floating-point array
    case 'double'
      for idx = 1:P
        if( isa(jArray(idx) , 'double' ))
          mArray(idx) = jArray(idx);
        else
          mArray(idx) = jArray(idx).doubleValue();
        end
      end
      mArray = double(mArray);
      
    case 'bu.util.Time'
      
      mArray = datenum( bu.util.Time.timeToDatevec( jArray )' );

    otherwise
      warning('Unknown type cast requested! Attempting to cast as double...' ); %#ok<WNTAG>
      for idx = 1:P
        mArray(idx) = jArray(idx).doubleValue();
      end
      mArray = double(mArray);
  end

end


