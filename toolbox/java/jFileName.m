function path = jFileName( str , flag )
tmp = java.io.File( str );
if( exist('flag','var') )
  if( flag )
    path = char(tmp.getName());
  else
    path = tmp.getName();
  end
else
  path = char(tmp.getName());
end