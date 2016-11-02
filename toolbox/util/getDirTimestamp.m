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
function [ filemap keys values ] = getDirTimestamp( directory , prefix , extension , dateformat )
%% Validate input arguments

filemap = java.util.TreeMap();
keys = [];
values = [];

% Determine class of object
% If the file parameter is not a java file object, convert it to one
directory = java.io.File(char(directory));
% Check if this is a valid directory
if( ~directory.isDirectory() ), return; end

% Convert prefix to java string
prefix = java.lang.String(prefix);
n_prefix = prefix.length();

% Check for a dot in the extension
extension = java.lang.String( extension );
if( extension.indexOf( '.' ) ~= 0 )
  extension = java.lang.String( ['.' char(extension)] );
end

% Get the length of the date format string
dateformat = java.lang.String( dateformat );
n_dateformat = dateformat.length();

%% Generate file listing
files = directory.listFiles();

% A treemap to store the files
filemap = java.util.TreeMap();

% Loop over files and find those that match the prefix and have the appropriate
% date formatting
for idx = 1:files.length
	
	time = jGetFileTimestamp( files(idx) );
	filemap.put( time , files(idx) );
	
%   
%   % Check for a match on the prefix
%   if( files(idx).getName().indexOf(prefix) == 0 )
%     
%     % Form suffix by eliminating prefix
%     suffix = files(idx).getName().substring(n_prefix);
%     
%     % We require an underscore separating the prefix and the date
%     if( suffix.charAt(0) ~= '_' ), continue; end
%     
%     % Take off the underscore
%     suffix = suffix.substring(1);
%     
%     % Take off the .mat extension
%     extIndex = suffix.indexOf( extension );
%     suffix = suffix.substring(0,extIndex);
%     
%     % Verify the length is the same as the date format string
%     if( n_dateformat ~= suffix.length() ), continue; end
%     
%     % Get time from format string
%     time = bu.util.Time.getTime( dateformat , suffix );
%     
%     % Add to filemap
%     filemap.put( time , files(idx) );
%     
%   end
end

if( nargout > 2 )
	keys = filemap.keySet().toArray( javaArray( 'bu.util.Time' , filemap.size() ) );
end
if( nargout == 3 )
	values = filemap.values().toArray( javaArray( 'java.io.File' , filemap.size() ) );
end

