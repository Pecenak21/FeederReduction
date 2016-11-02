%% Matlab Library
% Title: show
% Author: Bryan Urquhart
% Description:
%   Sets the visibility of the input argument to on. This is a shortcut to the
%   command set( h , 'visible' , 'on' );
%
function show( h )
  hl = get(0,'children');
  if( ismember(h,hl,'rows') )
    set(h,'visible','on');
  end
end