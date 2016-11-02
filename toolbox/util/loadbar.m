%% Matlab Library
%  Bryan Urquhart
%
%  Simple loading bar function to reduce code clutter
%
function loadbar( msg ,  numElem , nDots )

disp(msg);
java.lang.System.out.print('|');
for index = 1:numElem
  if( mod(index,floor(numElem/nDots)) == 0 ), java.lang.System.out.print('.'); end
end
java.lang.System.out.println('|');java.lang.System.out.print('|');