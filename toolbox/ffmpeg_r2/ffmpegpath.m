function ffmpegexe = ffmpegpath
%FFMPEGPAT   Returns the FFMPEG exe path

% Copyright 2013 Takeshi Ikuma
% History:
% rev. - : (05-15-2013) original release


if ispref('ffmpeg','exepath')
   ffmpegexe = getpref('ffmpeg','exepath');
else
   error('FFMPEG path not set. Run ffmpegsetup first. Download ffmpeg.exe from the internet.');
end
