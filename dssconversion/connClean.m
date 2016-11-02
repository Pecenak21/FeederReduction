function [conn, grounded] = connClean(str)
% clean up connection for OpenDSS

% clean input
str(str == ' ') = '';

if(any(strcmpi(str,{'wye','ln','yg','y','yng'}))), conn='wye';
elseif(any(strcmpi(str,{'ll','delta','d'}))), conn='delta';
else error('invalid connection type');
end

% Check if grounded
if strcmp(conn,'wye') % only consider grounded when wye conn is used.
	grounded = any(lower(str)=='g');
else
	grounded = 0;
end

end