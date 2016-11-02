%% Matlab Library
%  Bryan Urquhart
%
%  Description:
%    Merge two java arrays. Should be updated to handle n arrays...
%
function c = merge( a , b )
%% Check input arguments
if( class(a) ~= class(b) )
  error( 'Array type must be the same!' );
end

%% Perform Merge

% Allocate new array
c = javaArray( char(class(a(1))) , length(a) + length(b) );

% Merge arrays
index = 0;
for idx = 1:length(a)
  % Increment merged array index
  index = index + 1;
  % Add element to merged array
  c(index) = a(idx);
end

for idx = 1:length(b)
  % Increment merged array index
  index = index + 1;
  % Add element to merged array
  c(index) = b(idx);
end
