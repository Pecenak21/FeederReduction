%% Matlab Library
% Title: hide
% Author: Bryan Urquhart
% Description:
%   Sets the visibility of the input argument to off. This is a shortcut to the
%   command set( h , 'visible' , 'off' );
%
function hide( h )
  hl = get(0,'children');
  if( ismember(h,hl,'rows') )
    set(h,'visible','off');
  end
end