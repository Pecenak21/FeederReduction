%% Matlab Library
%
%  Title: Open matlab pool
%
%  Author: Bryan Urquhart
%
%  Description:
%    Collects some statements to open a matlab pool, close one if it exists.
%    This behavior is not complete, but is just collecting this code in one
%    spot so it doesn't take so many lines in other files that need to use it.
%
function lmMatlabpool( workers , quiet )

if( nargin == 0 )
  workers = 2;
end
if( nargin < 2 )
  quiet = 'quiet';
end

timer_ = tic;
isOpen = matlabpool('size') == workers;
if( ~isOpen )
  if( ~strcmpi('quiet',quiet) ), disp( 'Setting up matlab pool' ); end
  try matlabpool close force; catch e, end %#ok<NASGU>
  matlabpool(workers);
  if( ~strcmpi('quiet',quiet) ), fprintf( 'Pool setup complete. (%.2f sec.)\n',toc(timer_)); end
end