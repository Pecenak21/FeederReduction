function ffmpegtranscode(infile,outfile,varargin)
%FFMPEGTRANSCODE   Transcode multimedia file using FFmpeg
%   FFMPEGTRANSCODE(INFILE,OUTFILE) transcodes the input file specified by
%   the string INFILE using the H.264 video and AAC audio formats. The
%   transcoded data and outputs to the file specified by the string
%   OUTFILE. INFILE must be a FFmpeg supported multimedia file extension
%   (e.g., AVI, MP4, MP3, etc.) while the extension of OUTFILE is expected
%   to be MP4 (although it may output in other formats as well).
%
%   FFMPEGTRANSCODE(INFILE,OUTFILE,'OptionName1',OptionValue1,'OptionName2',OptionValue2,...)
%   may be used to customize the FFmpeg configuration:
%
%      Name    Description
%      ====================================================================
%      Range            Scalar or 2-element vector.
%                       Specifies the segment of INFILE to be transcoded.
%                       If scalar, Range defines the total duration to be
%                       transcoded. If vector, it specifies the starting
%                       and ending times.
%      Units            [{'seconds'}|'frames'|'samples']
%                       Time units for Range option. 'Frames' sets the time
%                       units to be the video frame index while 'samples'
%                       sets the units to be the audio sample index. Both
%                       'frames' and 'samples' options are defined with
%                       1-based indices.
%      FastSearch       ['on'|{'off'}]
%                       Enable or disable fast transcoding starting point
%                       search. With this option 'on', transcoding may
%                       begin substantially faster if starting time in
%                       Range option is substantially larger; however, it
%                       may result in less accurate starting position.
%      AudioCodec       [none|{copy}|wav|mp3|aac]
%                       Audio codec. If 'none', audio data would not be
%                       transcoded.
%      AudioSampleRate  Positive scalar
%                       Output audio sampling rate in samples/second.
%                       Only specify if needed to be changed.
%      Mp3Quality       Integer scalar between 0 and 9 {[]}
%                       MP3 encoder quality setting. Lower the higher
%                       quality. Empty uses the FFmpeg default.
%      AacBitRate       Integer scalar.
%                       AAC encoder's target bit rate in b/s. Suggested to
%                       use 64000 b/s per channel.
%      VideoCodec       [none|copy|raw|mpeg4|{x264}]
%                       Video codec. If 'none', video data would not be
%                       transcoded.
%      OutputFrameRate  Positive scalar
%                       Output video frame rate in frames/second.
%      InputFrameRate   Positive scalar
%                       Input video frame rate in frames/second. Altering
%                       the input frame rate effectively slows down or
%                       speeds up the video. This option is only valid for
%                       raw video format. Note that when both
%                       InputFrameRate and Range (with Units='seconds') are
%                       specified, Range is defined in the original frame
%                       rate.
%      VideoScale       Positive integer scalar or two-element vector
%                       Video size scaling factor. If scalar, the size of
%                       the output video is increased by the specified
%                       factor. If vector, it specifies the scaling factor
%                       as a ratio [num den]: num/den > 0 enlarges while
%                       num/den<0 shrinks the video frame size.
%      VideoCrop        4-element integer vector [left top right bottom]
%                       Video frame cropping/padding. If positive, the
%                       video frame is cropped from the respective edge. If
%                       negative, the video frame is padded on the
%                       respective edge.
%      VideoFillColor   ColorSpec
%                       Filling color for padded area.
%      VideoFlip        [horizontal|vertical|both]
%                       Flip the video frames horizontally, vertically, or
%                       both.
%      PixelFormat      One of format string returned by FFMPEGPIXFMTS
%                       Pixel format. Default to 'yuv420p' for Apple
%                       QuickTime compatibility.
%      x264Preset       [ultrafast|superfast|veryfast|faster|fast|medium|slow|slower|veryslow|placebo]
%                       x264 video codec options to trade off compression
%                       efficiency against encoding speed.
%      x264Tune         film|animation|grain|stillimage|psnr|ssim|fastdecode|zerolatency
%                       x264 video codec options to further optimize for
%                       input content.
%      x264Crf          Integer scaler between 1 and 51 {18}
%                       x264 video codec constant rate factor. Lower the
%                       higher quality, and 18 is considered perceptually
%                       indistinguishable to lossless. Change by ±6 roughly
%                       doubles/halves the file size.
%      Mpeg4Quality     Integer scalar between 1 and 31 {1}
%                       Mpeg4 video codec quality scale. Lower the higher
%                       quality.
%      DeleteSource     ['on'|{'off'}]
%                       Commands to delete all the input files at the
%                       completion.
%      ProgressFcn      ['none|{'default')|function handle]
%                       Callback function to display transcoding progress.
%                       For a custom callback, provide a function handle
%                       with form: progress_fcn(progfile,Nframes), where
%                       'progfile' is the location of the FFmpeg generated
%                       text file containing the transcoding progress and
%                       Nframes is the expected number of video frames in
%                       the output. Note that FFmpeg appends the new
%                       updates to 'progfile'. If set 'default', the
%                       transcoding progress is shown with a waitbar if
%                       video transcoding and no action for audio
%                       transcoding.
%
%   Example: Animation movie from a sequence of MATLAB plots:
%
%      % Generate one sinusoidal cycle with varying phase
%      t = linspace(0,1,1001);
%      phi = linspace(0,2*pi,21);
%      figure;
%      for n = 1:numel(phi)
%         plot(t,sin(2*pi*t+phi(n)))
%         print('-dpng',sprintf('test%02d.png',n)); % create an intermediate PNG file
%      end
%
%      % Create the MP4 file from the PNG files, animated at 5 fps
%      ffmpegtranscode('test%02d.png','sinedemo.mp4','InputFrameRate',5,...
%         'x264Tune','animation','DeleteSource','on');
%
%   References:
%      FFmpeg Home
%         http://ffmpeg.org
%      FFmpeg Documentation
%         http://ffmpeg.org/ffmpeg.html
%      FFmpeg Wiki Home
%         http://ffmpeg.org/trac/ffmpeg/wiki
%      Encoding VBR (Variable Bit Rate) mp3 audio
%         http://ffmpeg.org/trac/ffmpeg/wiki/Encoding%20VBR%20%28Variable%20Bit%20Rate%29%20mp3%20audio\
%      FFmpeg and AAC Encoding Guide
%         http://ffmpeg.org/trac/ffmpeg/wiki/AACEncodingGuide
%      FFmpeg and x264 Encoding Guide
%         http://ffmpeg.org/trac/ffmpeg/wiki/x264EncodingGuide
%      Xvid/Divx Encoding Guide
%         http://ffmpeg.org/trac/ffmpeg/wiki/How%20to%20encode%20Xvid%20/%20DivX%20video%20with%20ffmpeg
%      MeWiki X264 Settings
%         http://mewiki.project357.com/wiki/X264_Settings
%
%   See Also: FFMPEGSETUP, FFMPEGTRANSCODE

% Copyright 2013 Takeshi Ikuma
% History:
% rev. - : (06-19-2013) original release
% rev. 1 : (10-23-2013) Added 'VideoFlip' option
%                       Changed default x264Crf to 18

narginchk(2,inf);

% default options
opts = struct(...
   'Range',[],'Units','seconds',...
   'FastSearch','off',... % 'on'|{'off'}
   'AudioCodec','copy',... % none|copy|wav(pcm_s16le)|mp3|aac
   'AudioSampleRate',[],...
   'Mp3Quality',[],...
   'AacBitRate',[],...
   'VideoCodec','x264',... % none|copy|raw(AVI)|mpeg4|x264
   'OutputFrameRate',[],'InputFrameRate',[],...
   'VideoScale',[],...
   'VideoCrop',[],...
   'VideoFlip',[],...
   'VideoFillColor',[],...
   'PixelFormat','yuv420p',...
   'x264Preset','',... %ultrafast|superfast|veryfast|faster|fast|medium|slow|slower|veryslow|placebo
   'x264Tune','',... %film|animation|grain|stillimage|psnr|ssim|fastdecode|zerolatency
   'x264Crf',18,... % perceptually indistinguishable
   'Mpeg4Quality',1,... % best
   'DeleteSource','off','ProgressFcn','default');

% ‘-progress url (global)’
% ‘-r[:stream_specifier] fps (input/output,per-stream)’ Set frame rate
% ‘-ss position (input/output)’  When used as an input option (before -i), seeks in this input file to position. When used as an output option (before an output filename), decodes but discards input until the timestamps reach position. This is slower, but more accurate.
% ‘-t duration (output)’    Stop writing the output after its duration reaches duration. duration may be a number in seconds, or in hh:mm:ss[.xxx] form.

if ischar(infile)
   infile = cellstr(infile);
end

% check to make sure the input files exist
file = cellfun(@(f)which(f),infile,'UniformOutput',false);
I = cellfun(@isempty,file);
% if missing file contains '%' or if file can be located locally, let it pass ('which' function cannot resolve all the files)
if any(I) && any(cellfun(@(f)isempty(dir(f)),infile(I)) & cellfun(@(s)isempty(strfind(s,'%')),infile(I))) % ok 
   error('At least one of the specified files do not exist.');
else 
   infile(~I) = file(~I);
end

if ~(ischar(outfile) && size(outfile,1)==1)
   error('OUTFILE must be given as a string of characters.');
end

% if output file already exists, make sure it's not one of the input files
if exist(outfile,'file') 
   if any(strcmpi(infile,which(outfile)))
      error('INFILE and OUTFILE must be different.');
   end
end

% parse options (if given)
if nargin>2
   % gotta have even number of arguments
   if mod(nargin,2)>0
      error('Option parameters must be given as name-value pairs.');
   end
   
   % make sure names are valid
   names = varargin(1:2:end);
   if ~iscellstr(names)
      error('Option name must be given as a string of characters.');
   end
   
   % only process if value is non-empty
   vals = varargin(2:2:end);
   I = cellfun(@isempty,vals);
   names(I) = [];
   vals(I) = [];
   
   % determine which options are specified
   fnames = fieldnames(opts);
   [~,I] = ismember(lower(names),lower(fnames));
   if any(I==0)
      error('Invalid Option name.');
   end
   
   % set the options
   for n = 1:numel(I)
      opts.(fnames{I(n)}) = vals{n};
   end
end

% start constructing the command line
cmd = sprintf('%s -y', ffmpegpath());
%‘-y (global)’ Overwrite output files without asking.

% get progress file location
progfile = '';
if ~(isempty(opts.ProgressFcn) || strcmpi(opts.ProgressFcn,'none'))
   progfile = fullfile(tempdir,'ffmpegprogress.txt');
%    progfile = fullfile(pwd,'ffmpeg_progress.txt');
   if exist(progfile,'file')
      delete(progfile);
   end
   
   cmd = sprintf('%s -progress "%s"',cmd,progfile);
end

% Alter input frame rate
fs = opts.InputFrameRate;
if ~isempty(fs)
   if numel(infile)>1
      error('InputFrameRate can only be specified if there is only 1 input file.');
   end
   if ~(isnumeric(fs) && any(numel(fs)==[1 2]) && all(fs>0 & floor(fs)==fs & ~isinf(fs)))
      error('InputFrameRate must be given as a positive integer value or a pair of positive integers to represent a fraction, [num den].');
   end
   if numel(fs)==1
      cmd = sprintf('%s -r %d',cmd,fs);
   else
      cmd = sprintf('%s -r %d/%d',cmd,fs(1),fs(2));
   end
end

% specify input starting time if fastsearch is on
if strcmpi(opts.FastSearch,'on')
   if numel(infile)>1
      error('FastSearch can only be turned ''on'' if there is only 1 input file.');
   end
   [cmd,range,fs] = getrange(cmd,opts.Range,opts.Units,infile,opts.InputFrameRate,true);
   if range(1)>0
      opts.Units = 'seconds';
      opts.Range = range-range(1);
   end
end

% specify input file
for n = 1:numel(infile)
   cmd = sprintf('%s -i "%s"',cmd,infile{n});
end

% specify video codec
moreopt = true;
switch opts.VideoCodec
   case 'none'
      cmd = sprintf('%s -vn',cmd);
      moreopt = false;
   case 'copy'
      cmd = sprintf('%s -c:v copy',cmd);
      moreopt = false;
   case 'raw'
      cmd = sprintf('%s -c:v rawvideo',cmd);
   case 'mpeg4'
      q = opts.Mpeg4Quality;
      if ~(isnumeric(q) && numel(q)==1 && q>0 && q<32 && q==floor(q))
         error('Mpeg4Quality must be given as an integer between 1 and 31.');
      end
      
      % to-do: 2-pass (windows use NUL)
      % ffmpeg -y -i input.avi -c:v mpeg4 -vtag xvid -b:v 555k -pass 1 -an -f avi /dev/null
      % ffmpeg -i input.avi -c:v mpeg4 -vtag xvid -b:v 555k -pass 2 -c:a libmp3lame -b:a 128k output.avi
      cmd = sprintf('%s -c:v mpeg4 -q:v %d',cmd,q);
   case 'x264'
      q = opts.x264Crf;
      if ~(isnumeric(q) && numel(q)==1 && q>0 && q<52 && q==floor(q))
         error('x264Crf must be given as an integer between 1 and 51.');
      end
      
      cmd = sprintf('%s -c:v libx264 -crf %d',cmd,q);
      
      s = opts.x264Preset;
      if ~isempty(s)
         if ~(ischar(s) && size(s,1)==1)
            error('x264Preset must be given as a string of characters.');
         end
         cmd = sprintf('%s -preset %s',cmd,s);
      end
      
      s = opts.x264Tune;
      if ~isempty(s)
         if ~(ischar(s) && size(s,1)==1)
            error('x264Tune must be given as a string of characters.');
         end
         cmd = sprintf('%s -tune %s',cmd,s);
      end
      
      % to-do: 2-pass (windows use NUL)
      % ffmpeg -y -i input -c:v libx264 -preset medium -b:v 555k -pass 1 -an -f mp4 /dev/null && \
      % ffmpeg -i input -c:v libx264 -preset medium -b:v 555k -pass 2 -c:a libfdk_aac -b:a 128k output.mp4
      
   otherwise
      error('Unsupported VideoCodec specified.');
end

% If re-encoding video, add more options
if moreopt
   fs = opts.OutputFrameRate;
   if ~isempty(fs)
      if ~(isnumeric(fs) && any(numel(fs)==[1 2]) && all(fs>0 & floor(fs)==fs & ~isinf(fs)))
         error('OutputFrameRate must be given as a positive integer value or a pair of positive integers to represent a fraction, [num den].');
      end
      if numel(fs)==1
         cmd = sprintf('%s -r %d',cmd,fs);
      else
         cmd = sprintf('%s -r %d/%d',cmd,fs(1),fs(2));
      end
   end

   pix_fmt = opts.PixelFormat;
   if ~isempty(pix_fmt)
      if ~(ischar(pix_fmt) && size(pix_fmt,1))
         error('PixelFormat must be specified as a string of characters.');
      end
      cmd = sprintf('%s -pix_fmt %s',cmd,pix_fmt);
   end
   
   % set crop & scale filters if specified
   fcmd = setvideofilter(opts.VideoCrop,opts.VideoScale,opts.VideoFillColor,opts.VideoFlip);
   if ~isempty(fcmd)
      cmd = sprintf('%s -vf "%s"',cmd,fcmd);
   end
   
end

% specify audio codec
moreopt = true;
switch opts.AudioCodec
   case 'none'
      cmd = sprintf('%s -an',cmd);
      moreopt = false;
   case 'copy'
      cmd = sprintf('%s -c:a copy',cmd);
      moreopt = false;
   case 'wav'
      cmd = sprintf('%s -c:a pcm_s16le',cmd);
   case 'mp3'
      cmd = sprintf('%s -c:a mp3',cmd);
      
      q = opts.Mp3Quality;
      if ~isempty(q)
         if ~(isnumeric(q) && numel(q)==1 && q>=0 && q<10 && q==floor(q))
            error('Mp3Quality must be given as an integer between 0 and 9.');
         end
         cmd = sprintf('%s -q:a %d',cmd,q);
      end
   case 'aac'
      cmd = sprintf('%s -c:a aac -strict -2',cmd);
      
      b = opts.AacBitRate;
      if ~isempty(b)
         if ~(isnumeric(b) && numel(b)==1 && b>0 && b==floor(b) && ~isinf(b))
            error('x264Crf must be given as a positive integer.');
         end
         cmd = sprintf('%s -b:a %d',cmd,b);
      end
   otherwise
      error('Unsupported VideoCodec specified.');
end

if moreopt
   fs = opts.AudioSampleRate;
   if ~isempty(fs)
      if ~(isnumeric(fs) && numel(fs)==1 && fs>0 && floor(fs)==fs && ~isinf(fs))
         error('AudioSampleRate must be given as a positive integer value.');
      end
      cmd = sprintf('%s -ar %d',cmd,fs);
   end
end

% Trim option
try
   [cmd,range,fs] = getrange(cmd,opts.Range,opts.Units,infile,opts.InputFrameRate,false);
catch ME
   ME.rethrow;
end

% append output file name
cmd = sprintf('%s "%s"', cmd, outfile);

% set timer for progress
tobj = starttimer(opts.ProgressFcn,progfile,infile,range,fs);
try
   [s,r] = system(cmd);
   if s~=0
      % if failed, report the error
      error(r);
   end
catch ME
   if ~isempty(tobj)
      stop(tobj);
      delete(tobj);
      if exist(progfile,'file')
         delete(progfile);
         if strcmp(opts.ProgressFcn,'default')
            progfcn_default('delete');
         end
      end
   end
   ME.rethrow;
end
if ~isempty(tobj)
   stop(tobj);
   delete(tobj);
   
   if exist(progfile,'file')
      delete(progfile);
   end
   if strcmp(opts.ProgressFcn,'default')
      progfcn_default('delete');
   end
end

% delete all infiles if requested
if strcmpi(opts.DeleteSource,'on')
   for n = 1:numel(infile)
      if infileok(n)
         delete(infile{n});
      else
         [p,f,e] = fileparts(infile{n});
         if isempty(p), p = '.'; end
         files = dir(fullfile(p,['*' e]));
         files = {files.name};
         
         I = cellfun(@(file)isempty(sscanf(file,[f e])),files);
         files(I) = [];
         for m = 1:numel(files)
            delete(fullfile(p,files{m}));
         end
      end
   end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tobj = starttimer(progfcn,progfile,infile,range,fs)

tobj = [];

% if progress file not specified, progress not to be displayed
if isempty(progfile), return; end

try
   if isempty(range) && isempty(fs) % neither are specified
      [fs,T] = getframerate(infile);
   elseif isempty(range) % & ~isempty(fs)
      [fs0,T] = getframerate(infile);
      T = T*fs0/fs;
   else
      if isempty(fs) % & ~isempty(range)
         fs = getframerate(infile);
      end
      T = diff(range);
   end
   N = round(T*fs); % # of frames to be encoded
catch % audio coding does not produce trackable data so no progress bar to be shown
   return;
end

tobj = timer('ExecutionMode','fixedRate','Period',1,'TasksToExecute',inf);
if strcmpi(progfcn,'default')
   set(tobj,'TimerFcn',@(~,~)progfcn_default(progfile,N));
   progfcn_default('setup');
elseif isa(progfcn,'function_handle')
   set(tobj,'TimerFcn',{progfcn,progfile,N});
else
   error('ProgressFcn must be given as a function handle object.');
end

start(tobj);

end

function progfcn_default(progfile,N)
persistent pos
persistent h

switch progfile
   case 'setup'
      pos = 0;
      h = waitbar(0,'Searching for the starting frame...','WindowStyle','modal',...
         'Name',mfilename,'CloseRequestFcn',{});
      drawnow;
   case 'delete'
      delete(h);
   otherwise
      
      fid = fopen(progfile,'r');
      if fid<0, return; end % file not ready
      txt = '';
      try
         fseek(fid,pos,-1); % go to the last position
         txt = fscanf(fid,'%c',inf); % read all the text to the end
         pos = ftell(fid); % save the last position
         fclose(fid);
      catch % just in case
         fclose(fid);
      end
      
      toks = regexp(txt,'frame=(\d+)\n','tokens');
      if ~isempty(toks)
         val = str2double(toks{end}{1});
         if val>0
            waitbar(val/N,h,'Video transcoding in progress...');
            drawnow;
         end
      end
end
end


function [cmd,range,fs_in] = getrange(cmd,range,units,infile,fs_in,input)

% if range not specified, nothing to do
if isempty(range), return; end

N = numel(range);
if ~(isnumeric(range) && any(N==[1 2]) && any(range>=0 & ~isnan(range) & ~isinf(range)))
   error('Range must be give as a positive scalar or 2-element vector.');
end
if N==2 && range(2)<=range(1)
   error('Range values must be increasing.');
end

if ~(ischar(units) && size(units,1)==1)
   error('Units must be given as a string of characters.');
end

if strcmpi(units,'seconds') % adjust if input frame rate changed
   if N==1 && range(1)==0
      error('Span of Range must be greater than 0.');
   end
   if ~isempty(fs_in) % if input rate deliverately changed, modify the range
      fs0 = getframerate(infile);
      range = range*fs0/fs_in;
   end
   if N==1
      range = [0 range];
   end
else
   if N==1 && range(1)==1
      error('Span of Range must be greater than 0.');
   elseif N==2 && range(1)==0
      error('Sample/frame based Range must be given with 1-based indices.');
   end
   
   if strcmpi(units,'frames') % video frame
      if isempty(fs_in)
         fs_in = getframerate(infile);
      end
   elseif strcmpi(units,'samples')
      fs_in = getsamplerate(infile);
   else
      error('Invalid Units value specified.');
   end
   
   if N==1
      range = [0 range/fs_in];
   else
      range = (range(:).'-[1 0])/fs_in;
   end
end

if range(1)>0
   cmd = sprintf('%s -ss %0.6f',cmd,range(1));
end
if ~input
   cmd = sprintf('%s -t %0.6f',cmd,diff(range));
end

end

function [fs0,T] = getframerate(infile)

info = ffmpeginfo(infile);
found = false;
for m = 1:numel(info)
   s = info(m).streams;
   for n = 1:numel(s)
      found = strcmp(s(n).type,'video');
      if found
         break;
      end
   end
end
if ~found
   error('Input file does not contain any video stream.');
end

c = s(n).codec;
if isempty(c.fps)
   if isempty(c.tbr)
      if isempty(c.tbn)
         if isempty(c.tbc)
            error('Input file does not report its frame rate.');
         else
            fs0 = c.tbc;
         end
      else
         fs0 = c.tbn;
      end
   else
      fs0 = c.tbr;
   end
else
   fs0 = c.fps;
end

if nargout>1
   T = info(m).duration;
end

end

function [fs0,T] = getsamplerate(infile)

info = ffmpeginfo(infile);
found = false;
for m = 1:numel(info)
   s = info(m).streams;
   for n = 1:numel(s)
      found = strcmp(s(n).type,'audio');
      if found
         break;
      end
   end
end
if ~found
   error('Input file does not contain any audio stream.');
end

fs0 = s(n).codec.samplerate;
if isempty(fs0)
   error('Input file''s audio stream does not specify sampling rate.');
end

if nargout>1
   T = info(m).duration;
end

end

function fcmd = setvideofilter(crop,scale,color,flip)

fcmd = '';
N = numel(crop);
if N>0
   if ~(isnumeric(crop) && numel(crop)==4 && all(crop==floor(crop) & ~isinf(crop)))
      error('VideoCropping must be specified as [top left bottom right] in pixels.');
   end

   pad = max(-crop,0);
   crop = max(crop,0);
   
   if any(crop>0)
      wadj = sum(crop([1 3]));
      hadj = sum(crop([2 4]));
      fcmd = sprintf('%scrop=in_w-%d:in_h-%d:%d:%d',fcmd,wadj,hadj,crop(1),crop(2));
   end
   if any(pad>0)
      if ~isempty(fcmd)
         fcmd = sprintf('%s, ',fcmd);
      end

      wadj = sum(pad([1 3]));
      hadj = sum(pad([2 4]));
      fcmd = sprintf('%spad=in_w+%d:in_h+%d:%d:%d',fcmd,wadj,hadj,pad(1),pad(2));
      
      N = numel(color);
      if N>0
         if ischar(color) 
            if N>1
               val = find(strcmpi(color,{'black','red','green','yellow','blue','magenta','cyan','white'}));
            else
               val = find(color=='krgybmcw',1);
            end
            if isempty(val)
               error('Invalid color specifier for VideoFillColor option value.');
            end
            color = bitget(val-1,1:3);
         elseif ~(isnumeric(color) && N==3 && all(color>=0 & color<=1))
            error('Invalid color specifier for VideoFillColor option value.');
         end
         color = dec2hex(floor(color*255),2);
         fcmd = sprintf('%s:color=%s',fcmd,reshape(color.',1,6));
      end
   end
end

N = numel(scale);
if ~isempty(scale)
   if ~(isnumeric(scale) && any(N==[1 2]) && all(scale>0&~isinf(scale)&scale==floor(scale)))
      error('VideoScaling must be specified as a positive integer value or a pair of positive integers to represent a fraction [num den]');
   end
   
   if ~isempty(fcmd)
      fcmd = sprintf('%s, ',fcmd);
   end
   
   if N==1
      fcmd = sprintf('%sscale=w=%d*iw:h=-1',fcmd);
   else %if N==2
      fcmd = sprintf('%sscale=w=%d/%d*iw:h=-1',fcmd,scale(1),scale(2));
   end
end

if ~isempty(flip)
   if ~(ischar(flip) && isrow(flip))
      error('VideoFlip must be specified as a string.');
   end
   I = find(strcmpi(flip,{'horizontal','vertical','both'}));
   if isempty(I)
      error('Invalid VideoFlip option given.');
   end
   if mod(I,2)
      if ~isempty(fcmd)
         fcmd = sprintf('%s, ',fcmd);
      end
      fcmd = sprintf('%shflip',fcmd);
   end
   if I>1
      if ~isempty(fcmd)
         fcmd = sprintf('%s, ',fcmd);
      end
      fcmd = sprintf('%svflip',fcmd);
   end   
end

end
