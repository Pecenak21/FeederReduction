function path = jFileStr( str )
tmp = java.io.File( str );
path = char(tmp.getPath());