%% Matlab Library
%  Bryan Urquhart

% This function performs inverse distance weighting given a point to
% interpolate, a set of points to interpolate from, and the order of the
% distance weighting
%
% Input:
% x - point to be interpolated
% z - interpolation set, column vector with third entry as magnitude
% a - order of distance weighting
%
%
function u = idistw( x , z , a )
%% Input checking

% Verify the size of x
if( size(x,1) ~= 1 || size(x,2) ~= 2 )
  error( 'The input parameter x must be a 1 x 2 input' );
end

% Verify the size of y
if( size(z,2) ~= 3 )
  error( 'The input parameter y must be a n x 3 input' );
end

%% Distance computation
%  This section computes the distances from each point in z to the point x

% Get the interpolation set size
N = size(z,1);

% Allocate a distance vector
distance = zeros(N,1);

% Loop over interpolation set
for i = 1:N
  % Compute the distance as the 2 or euclidean norm
  distance(i) = sqrt( (z(i,1)-x(1,1))^2 + (z(i,2)-x(1,2))^2 );
end

%% Compute the weight normalization

% Initialize weight normalization
k = 0;

% Loop over interpolation set
for i = 1:N
  k = k + 1 / distance(i)^a;
end

%% Compute the inverse distance weighting

% Initialize IDW result to zero
u = 0;

% Loop over interpolation set
for i = 1:N
  u = u + z(i,3) / distance(i)^a;
end

u = u / k;
