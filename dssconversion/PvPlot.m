function h = PvPlot(circuit, varargin)

%% Create the figure and add items to it
% figure and it's position
% sz = [820 530];
h = figure('Name',['CircuitVisualizer: ' circuit.circuit.Name],'toolbar','figure','numbertitle','off');
p = get(h,'pos');
% p(2) = p(2)-(sz(2)-p(4));
% p(3:4) = sz;
set(h,'pos',p);

% axes on the left half
cs.axes = axes('XTick',[],'YTick',[]);
circuit.n = setdiff(fieldnames(circuit),{'buslist','basevoltages','n'});
if(isfield(circuit,'buslist'))
	circuit.n = ['bus' circuit.n(:)'];
end
circuit.cs = cs;

%% perform data calculations

% build a struct lookup table for node coords
bl = circuit.buslist;
bl.id = fnSanitize(bl.id);
for i_=1:length(bl.id)
	bn = ['b' bl.id{i_}];
	l.(bn).xy = bl.coord(i_,:);
	l.(bn).ind = i_;
	l.(bn).dev = {};
end

% build a hash lookup table for sections
sl = circuit.line;
buslist = strcat('b',stripPhases({sl.bus1; sl.bus2}'));
mask = any(strcmp(buslist,'b'),2);
buslist = mat2cell(buslist,ones(1,size(buslist,1)),size(buslist,2));
buslist(mask) = [];

namelist = strcat('line_',sl(~mask).Name);
% and, if applicatble, transformers
if(isfield(circuit,'transformer'))
	tl = circuit.transformer;
	namelist = vertcat(namelist,strcat('transformer_',{tl.Name}'));
	nbl = cell(length(tl),1);
	for i_=1:length(tl)
        nbl{i_} = stripPhases(tl(i_).Bus);
	end
	nbl = cellfun(@(x)strcat('b',x),nbl,'UniformOutput',false);
	buslist = [buslist; nbl];
end
% actually add the lookup
buslist = cellfun(@fnSanitize,buslist,'UniformOutput',false);
for i_=1:length(buslist)
	bl = buslist{i_};
	bn = namelist{i_};
	sec.(bn).bus = bl;
	sec.(bn).ind = i_;
	l.(bl{1}).dev{end+1} = bn;
	l.(bl{2}).dev{end+1} = bn;
end

% all the others
o_name = setdiff(circuit.n,{'buslist','bus','line','transformer'});
pmap = struct('load','bus1','generator','bus1','circuit','bus1','capacitor','Bus1','regcontrol','bus');
missing_items = 0;
for objtype = o_name; objtype = objtype{1};
	ol = circuit.(objtype);
	if(~isfield(pmap,objtype)), continue; end;
	namelist = strcat([objtype,'_'],{ol.Name});
	bus_fn = pmap.(objtype);
	bl = strcat('b',stripPhases(ol.(bus_fn)));
	if(~iscell(bl)), bl = {bl}; end
	for i_ = 1:length(ol)
		try % only mark buses that we know about
			l.(bl{i_}).dev{end+1} = namelist{i_};
		catch
			missing_items = missing_items +1;
		end
	end
end

% save the data for later use
circuit.bl = l;
circuit.sl = sec;

%% draw the whole circuit
% call the addToPlot subfunction to do the work
circuit.plots{1} = addToPlot(circuit, false, fieldnames(sec));

%% add loads
bl = circuit.buslist;
load = [circuit.load];
ld = struct();
for i = 1:length(load)
    ld(i).kVA = sqrt(load(i).Kw^2+load(i).Kvar^2);
    ld(i).bus1 = load(i).bus1;
    ld(i).Name = load(i).Name;
end
% add useful fields
ld().lon = [];
ld().lat = [];
% update lat lon
for i = 1:length(ld)
	% check for .1.2.3 and get rid of it from the bus name
	z = ld(i).bus1 == '.'; 
	if ~isempty(find(z,1))
		z = find(z,1,'first');
		z = ld(i).bus1(1:z-1);
	else
		z = ld(i).bus1;
	end
	[x, y] = ismember(z,bl.id);
	ld(i).lon = bl.coord(y,1);
	ld(i).lat = bl.coord(y,2);
end
cmap = colormap(jet(256));
for i = 1:length(ld)
    if ceil( log10(ld(i).kVA)/log10(1000)*255 ) < 0
        disp('Load too small to be plotted');%assume smallest possible value
        plot([ld(i).lon]',[ld(i).lat]','.r','markersize',log(1.001)*15,'linewidth',2,...
			'color',cmap(ceil( 1.001 ),:) );
        continue
    elseif ceil( log10(ld(i).kVA)/log10(1000)*255 ) > 256
        disp('Load bigger than color array');%assume biggest possible value
		plot([ld(i).lon]',[ld(i).lat]','.r','markersize',log(ld(i).kVA)*15,'linewidth',2,...
			'color',cmap(ceil( 256 ),:) );
        continue
    end
		plot([ld(i).lon]',[ld(i).lat]','.r','markersize',log(ld(i).kVA)*15,'linewidth',2,...
			'color',cmap(ceil( log10(ld(i).kVA)/log10(1000)*255 ),:) );
end
% xlabel('Longitude','fontsize',20);
% ylabel('Latitude','fontsize',20);
set(gca,'fontsize',20);
set(gca,'ydir','normal');
% Rescale color bar
d = log10([ld.kVA]);
mn = min(d(:));
rng = max(d(:))-mn;
d = 1+63*(d-mn)/rng; % Self scale data
hC = colorbar;
L = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000];
% Choose appropriate
% or somehow auto generate colorbar labels
l = 1+255*(log10(L)-mn)/rng; % Tick mark positions
set(hC,'Ytick',l,'YTicklabel',L,'fontsize',20);
ylabel(hC,'Load rating, [kVA]','fontsize',20);

%% Add items to the Plot
	function lo = addToPlot(c,isNode,oname,clr,line_width)
		% addToPlot is responsible for looking up and plotting new data
		% it takes a circuit data structure (including control handles),
		% a flag indicating whether the objects are lines or nodes, the
		% object names, and an optional color to use.
		% Returns handles to the graph objects that were created
		
		% first get the axes and make sure we'll actually be able to add to
		% them.
		ah = c.cs.axes;
		set(ah,'Next','add');
		diag_l = hypot(diff(get(ah,'YLim')),diff(get(ah,'XLim')))/751; % 751 is approximately the number of pixels in the diagonal of the default (and therefore, probably smallest) size of the graph area, so this is a rough measure of units/pixel of the graph

		% default to black if no color was given
		if(nargin < 4 || isempty(clr)); clr = 'black'; end;
		% default to line width of 1:
		if(nargin < 5 || isempty(line_width)); line_width = 1; end;
		% make sure the object names are in the format we want
		if(~iscell(oname)); oname = {oname}; end
		oname = fnSanitize(oname);
		
		% allocate an array to hold the handles for output.
		lo = zeros(length(oname),1);
		% lookup the data in the cached table and plot them one at a time
		if(isNode)
			for ii=1:length(oname)
				coords = c.bl.(oname{ii}).xy;
				lo(ii) = plot(ah,coords(1),coords(2),'x','Color',clr,'UserData',{isNode,oname{ii}},'LineWidth',line_width,'MarkerSize',max(6,line_width*6));
			end
		else
			for ii=1:length(oname)
				bn = c.sl.(oname{ii}).bus;
				coords = vertcat(l.(bn{1}).xy,l.(bn{2}).xy);
				seg_l = hypot(diff(coords(:,1)),diff(coords(:,2)))/diag_l; % this should be the pixel length of this segment
				lo(ii) = plot(ah,coords(:,1),coords(:,2),'Color',clr,'UserData',{isNode,oname{ii}},'LineWidth',line_width,'MarkerSize',max(6,line_width*6));
				if(seg_l < 4)
					set(lo(ii),'Marker','x');
				end
			end
		end
    end
end