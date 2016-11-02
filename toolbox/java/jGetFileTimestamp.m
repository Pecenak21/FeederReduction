%% Solar Resource Assessment
%  Matlab Library Development
%  Author: Bryan Urquhart
%  Date:   July 12, 2011
% 
%
%  Description:
%    This function obtains time stamped file names with a given core name
%    followed by an underscore and a time stamp
%
function time = jGetFileTimestamp( filename )
%% Process input arguments


%% Convert filename

file = java.io.File(char(filename));
name = char(file.getName());
index = regexp( name , '\d{14}' );
tStr = file.getName().substring(index-1,index+14-1);
time = bu.util.Time.YYYYMMDDHHMMSS(tStr);

