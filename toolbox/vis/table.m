%% Matlab Library
%
% Author: Bryan Urquhart
%
% Title: Generate Table
%
% Description:
%   Make a table with data
%
%
function h = table( data , rowheading , colheading , varargin )
%% Process Input arguments

% rowheading = {'\Deltat','Stein','Iman'};
% colheading = {};
% data = [interval' vi_stein' vi_iman'];

nrowheadings = size( rowheading , 1 );
ncolheadings = size( colheading , 2 );

[m n] = size(data);

% Default varargs
visible = 'on';

if( ~isempty( varargin ) )
	args = argHandler( varargin );
	for idx = 1:size(args,1)
		switch( args{idx,1} )
			case 'visible'
				visible = args{idx,2};
		end
	end
end
 
 %%
nrows = m+nrowheadings;
ncols = n+ncolheadings;
 
%if( ~isempty( rowheading ) ), nrows = nrows+1; end
%if( ~isempty( colheading ) ), ncols = ncols+1; end

 %% Set up figure
 rowsize = 55;
 colsize = 220;
 rowbuf = 200;
 colbuf = 200;
 
 h = figure;
 set( h , 'visible' , visible );
 set( h , 'color' , [1 1 1]);
 set( h , 'Position' , [300 300 colsize*ncols+colbuf rowsize*ncols+rowbuf] );
 set( h , 'InvertHardcopy' , 'off' );
 set( h , 'PaperPositionMode' , 'auto' );
  
 ha = axes;
 
%% Draw grid
 
for ridx = 0:nrows
	linewidth = 1;
	if( ridx == nrowheadings ), linewidth = 2.5; end
	if( ridx == 0 | ridx == nrows ),        linewidth = 4; end
	line( [0 ncols] , [ridx ridx] , 'color' , [0 0 0] , 'linewidth' , linewidth );
end

for cidx = 0:ncols
	if( cidx == ncolheadings ), linewidth = 2.5; end
	if( cidx == 0 | cidx == ncols ),        linewidth = 4; end
	line( [cidx cidx] , [0 nrows] , 'color' , [0 0 0]  , 'linewidth' , linewidth);
end

%% Populate table

for ridx = 1:nrows
	for cidx = 1:ncols
		
		[ridx cidx]
		if( ridx < nrowheadings+1)
			str = rowheading{ridx,cidx};
			text('Position',[cidx+1/8-1 ridx-1/2+nrowheadings-1],'String',str,'Color',[0 0 0],'FontSize',14,'FontWeight','bold');
			continue;
		end
		
		if( cidx < ncolheadings+1)
			str = colheading{ridx-nrowheadings,cidx};
			text('Position',[cidx+1/8-1 ridx-1/2+nrowheadings-1],'String',str,'Color',[0 0 0],'FontSize',14,'FontWeight','bold');
			continue;
		end
		
		str = sprintf( '%5g', data(ridx-nrowheadings,cidx-ncolheadings) );
		text('Position',[cidx+1/8-1 ridx-1/2+nrowheadings-1],'String',str,'Color',[0 0 0],'FontSize',14,'FontWeight','bold');
	end
end

 %% Set up axes
 
 set( ha , 'xlim' , [-1 ncols+1] );
 set( ha , 'ylim' , [-1 nrows+1] );
 set(ha,'ydir','reverse');
 axis off;
 box off;
 
 %set(ha, 'Position', get(ha, 'OuterPosition') - get(ha, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
