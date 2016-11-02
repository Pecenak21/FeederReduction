function t = nowUTC()
% like now() but in UTC instead of local clock time

% use the Java interface to get the UTC rather than local time.
% convert milliseconds to days, then add the offset for the Unix Epoch:
% 719529 = datenum('January 1, 1970');
t = java.lang.System.currentTimeMillis/1000/3600/24 + 719529;
