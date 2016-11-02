%% Matlab Visualization Library
%
%  Author: Bryan Urquhart
%  Date:   2012-12-14
%
%  Description:
%    Plots a surface using patch elements
%
%
%%
function hs = visSurfPatch( varargin )
%% Process Input Argument

switch( nargin )
  case 1
    z = varargin{1};
    [m n] = size(z);

    x = 1:n;
    y = 1:m;
    
    cmap = jet(256);
    
  case 3
    x = varargin{1}(1,:);
    y = varargin{2}(:,1);
    z = varargin{3};
    cmap = jet(256);
    
    [m n] = size(z);
    
    if( numel(x) ~= n || numel(y) ~= m )
      error( 'Input arguments have inconsistent dimensions');
    end
    
  case 4
    x = varargin{1}(1,:);
    y = varargin{2}(:,1);
    z = varargin{3};
    cmap = varargin{4};
    
    [m n] = size(z);
    
    if( numel(x) ~= n || numel(y) ~= m )
      error( 'Input arguments have inconsistent dimensions');
    end
    
    if( size(cmap,2) ~=3 )
      error( 'Colormap must have 3 columns' );
    end
    
  otherwise
    error( 'Unsupported number of arguments' );
end

q = size(cmap,1);
zlim = [min(z(:)) max(z(:))];

%% Make patches

flag3D = true;

XData = nan(3,2*m*n);
YData = nan(3,2*m*n);
ZData = nan(3,2*m*n);
CData = nan(3,2*m*n);
cnt = 0;

for i = 1:m-1
  for j = 1:n-1
    
    % Construct triangle 1
    xdata = [ x(j) x(j+1) x(j)   ; x(j) x(j+1) x(j+1) ]';
    ydata = [ y(i) y(i+1) y(i+1) ; y(i) y(i)   y(i+1) ]';
    zdata = [ z(i,j) z(i+1,j+1) z(i+1,j) ; z(i,j) z(i,j+1) z(i+1,j+1) ]';
   
    if(sum(isnan(zdata(:)))), continue; end
 
    cidx = getColorIndex(zdata,zlim,q);
    %hs = patch(xdata,ydata,zdata,cidx,'EdgeColor','none','FaceColor',cmap(cidx(1),:));
    
    cnt = cnt + 2;
    
    XData(:,cnt-1:cnt) = xdata;
    YData(:,cnt-1:cnt) = ydata;
    ZData(:,cnt-1:cnt) = zdata;
    CData(:,cnt-1:cnt) = cidx;
    
  end
end

hs = patch(XData,YData,ZData,CData,'EdgeColor','none');
set(gcf,'colormap',cmap);

%% Set the camera position

xrange = max(x(:)) - min(x(:));
yrange = max(y(:)) - min(y(:));
zrange = max(z(:)) - min(z(:));
cpos(1) = x(1)-xrange;
cpos(2) = y(1)-yrange;
cpos(3) = max(z(:)) + 2*zrange;

%set(gca,'CameraPosition',cpos);

end

%% Support functions
function idx = getColorIndex( val , zlim , n )
  idx = ceil( ( val - zlim(1) ) / ( zlim(2) - zlim(1) ) * n );
  idx = max(1,idx); idx = min(255,idx);
end