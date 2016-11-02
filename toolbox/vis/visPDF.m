%% Solar Resource Assessment
%  Copper Mountain Library
%
%  Title: Generate the kt PDF for a day
%
%  Author: Bryan Urquhart
%
%  Description:
%    Using the power data stored within the pwr variable along with the clear
%    sky structure, generate a PDF of power.
%
function [h ha] = visPDF( f , varargin )
%% Process Input Arguments

% csk = [];
% cskfile = [];
% dat = [];
% cfg = [];

visible = 'on';
zenithmax = 86;
nbins = 500;
fontsize = 14;
fontweight = 'normal';
linewidth = 1;

% Process varargin
if( ~isempty( varargin ) )
  args = argHandler( varargin );
  for idx = 1:size(args,1)
    switch( args{idx,1} )
      case 'visible'
        visible = args{idx,2};
      case 'zenithmax'
        zenithmax = args{idx,2};
      case 'nbins'
        nbins = args{idx,2};
			case 'fontsize'
				fontsize = args{idx,2};
			case 'fontweight'
				fontweight = args{idx,2};
			case 'linewidth'
				linewidth = args{idx,2};
			case 'color'
				color = args{idx,2};
    end
  end
end

%% Plot

% Open figure
h  = figure;
ha = axes;

% Set up figure
set(h,'Visible',visible);
set(h, 'color' , [1 1 1]);
set(h, 'paperpositionmode', 'auto');
set(h, 'inverthardcopy', 'off');

set(ha,'FontSize',fontsize);
set(ha,'fontweight',fontweight);

x = f.bins(:,2);

% Plot the pdf
line( x , f.binval/f.integral , 'linewidth' , linewidth , 'color' , color );
hold on;
%plot( ha , x , f.binval_filt/f.integral , '--' , 'linewidth' , 2 , 'color' , 'r' );

% Set up the axes
xmin = min( f.bins(:) );
xmax = max( f.bins(:) );
set( ha , 'Xlim' , [xmin xmax] );
% xlim = get( ha , 'Xlim' );
% xmin = xlim(1);
% xmax = xlim(2);
%xmin = 0;
%xmax = ceil( max( f.bins(:,3) ) * 10 ) / 10;
%set( ha , 'XTick' , xmin:0.2:xmax );

ylim = get( ha , 'Ylim' );
ymin = ylim(1);
ymax = ylim(2);
% ymin = 0;
% ymax = ceil( max( ktpdf.binval/ktpdf.integral ) * 10 ) / 10;
%set( ha , 'YTick' , ymin:0.2:ymax );

% Label the axes
xlabel( 'Abscissa [ ]','FontSize',fontsize,'fontweight',fontweight );
ylabel( 'Probability Density [ ]','FontSize',fontsize,'fontweight',fontweight );

% Legend
%legend('orginal','filtered','Location','Best');

% Text insertion
% Compute text position
x = xmin + (xmax-xmin)*.45;
y = ymax - (ymax-ymin)*0.05;
%text('Position',[x y],'String',['bin width: ' num2str(mean(f.binwidth))],'Color',[0.6 0.8 1]*0.4,'FontSize',10,'FontWeight','bold','FontName','Arial');
