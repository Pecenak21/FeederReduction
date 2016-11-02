%% Matlab Library
%  Bryan Urquhart

%  Path includes for sub folders used in analysis

function rm_asv( directory , varargin )
%% Process Input Arguments

visible = 'off';

if( ~isempty( varargin ) )
  args = argHandler(varargin);
  for idx = 1:size(args,1)
    switch( args{idx,1} )
      case 'visible'
        visible = args{idx,2};
    end
  end
end

%% Folders

dir = java.io.File( char(directory) );

sub_rm_asv( dir , visible );

end

function sub_rm_asv( dir , visible )
  listing = dir.listFiles();

  for idx = 1:listing.length
    if( listing(idx).isFile() )
      if( listing(idx).getName().endsWith('.asv') )
        if( strcmpi(visible,'on') )
          disp( [ 'Deleting ' char( listing(idx) ) ] );
        end
        listing(idx).delete();
      end
    end
    if( listing(idx).isDirectory() )
      sub_rm_asv( listing(idx) , visible );
    end
  end
end