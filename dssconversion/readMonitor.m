function m = readMonitor(d)
% Reads in the byte stream from an OpenDSS monitor object
% Refer to
% http://sourceforge.net/projects/electricdss/forums/forum/861976/topic/4765650
%% Read header
% Signature: 32-bit integer (should be 43756) 
m.sig = typecast(d(1:4),'int32');
if m.sig ~= 43756, 
    error('ByteStream did not contain expected signature'); 
end
% Version: 32-bit integer 
m.ver = typecast(d(5:8),'int32');
% Recordsize: 32-bit integer (bytes each record)
m.size = typecast(d(9:12),'int32');
% Mode: 32-bit integer (monitor mode)
m.mode = typecast(d(13:16),'int32');
% Header String: 256-byte ANSI characters (fixed length)
m.header = native2unicode(d(17:272));
%% Read Reacords
% Channels repeat every m.size + 2 times (+2 for hour and sec records)
out = typecast(d(273:end),'single');
out = reshape(out, m.size+2, [])';
m.data = out;
end