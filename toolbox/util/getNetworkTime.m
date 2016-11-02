%% Solar Resource Assessment
%  Image Acquisition
%
%  Author: Bryan Urquhart
%
%  Description:
%    Accesses NIST site and get the current network time
%
%
%% Check network time server
function time = getNetworkTime()

% Set up NTP info
port = 13;
machine = 'time.nist.gov';

time = -1;
try
  
  count = 0;
  
  while( true )
    
    % Set up the socket
    sock = java.net.Socket(machine,port);

    % Set the input stream
    istream = sock.getInputStream();

    % Set the input stream reader
    isreader = java.io.InputStreamReader( istream ) ;

    % Set up the buffered reader
    reader = java.io.BufferedReader( isreader );

    % Read a blank line, then the time string
    reader.readLine();
    timestr = reader.readLine();
    
    % Check for blank time returns
    if( timestr.length() < 20 )
      pause(1);
      % Increment counter
      count = count + 1;
      continue;
    end
    
    % Parse the time string for the time
    time = timestr.split('[ ]');
    yymmdd = char(time(2).split('[-]'));
    hhmmss = char(time(3).split('[:]'));
    
    y  = str2double( [ '20' yymmdd(1,:) ] );
    m  = str2double( yymmdd(2,:) );
    d  = str2double( yymmdd(3,:) );
    hh = str2double( hhmmss(1,:) );
    mm = str2double( hhmmss(2,:) );
    ss = str2double( hhmmss(3,:) );
    
    time = datenum( [ y m d hh mm ss ] );   
    
    % Check for blank time returns
    if( timestr.length() > 20 || count > 5 )
      break;
    end
    
  end
  
catch e %#ok<*NASGU>

end

try
  sock.close();
  istream.close();
  isreader.close();
  reader.close();
catch e
end