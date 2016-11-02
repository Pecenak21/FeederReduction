function [h, cs] = circuitVisualizer(circuit, varargin)
% circuitVisualizer is a simple graphical interface to help inspect which
% parts of a circuit connect to each other.
%
% To run it, you need a dss circuit structure (such as produced by the
% dssconversion function, or used by dsswrite).  Returns the figure handle
% of the resulting window, which you can just ignore for most purposes.

% to allow outside calling of our subfunctions, we support a syntax where
% if the first argument is a string, we try to call our subfunction with
% that name and the remaining arguments.  This is sometimes useful for
% callback functions although I think I haven't actually used it in this
% case.
if(ischar(circuit))
	eval([circuit '(varargin{:});']);
	return;
% elseif(isstruct(circuit) && isfield(circuit,'load'))
elseif(isstruct(circuit)) % js: some circuits do not have loads, so I deactivated that check
	%This is just a quick check that we can deal with the input data.  In
	%this case we let control fall through into the rest of the function.
else
	error('circuitVisualizer:invalidInput','Must input circuit struct');
end

% We can only draw circuits with coordinate data
if(~isfield(circuit,'buslist'))
	error('circuitVisualizer:noLocationData','Circuits must contain location data to be able to visualize');
end

%% Create the figure and add items to it
% figure and it's position
sz = [820 530];
h = figure('Name',['CircuitVisualizer: ' circuit.circuit.Name],'toolbar','figure','numbertitle','off'); 
p = get(h,'pos');
p(2) = p(2)-(sz(2)-p(4));
p(3:4) = sz;
set(h,'pos',p,'color','w');
% Normally I would go ahead and set the resize function here, but in newer
% versions of MATLAB, it causes problems when I go to get the java
% component for the edit field (the resizefcn gets called, and we haven't
% assigned the requisite data handle for it to use yet at that point), so
% I'm moving it down to be set after that.

% axes on the left half
cs.axes = axes('XTick',[],'YTick',[]);
set(cs.axes,'NextPlot','add','DataAspectRatio',[1 1 1],'box','on');
% menu of circuit types
cs.types = uicontrol('Style','listbox');
circuit.n = setdiff(fieldnames(circuit),{'buslist','basevoltages','n'});
if(isfield(circuit,'buslist'))
	circuit.n = ['bus' circuit.n(:)'];
end
set(cs.types,'String',circuit.n,'Value',1);
set(cs.types,'Callback',@(o,e)typeCB(o));
% menu of object names
cs.objects = uicontrol('Style','listbox');
set(cs.objects,'Callback',@(o,e)typeCB(o));
% field with object name
cs.name = uicontrol('Style','edit');
% to do autocomplete well, we would want to add text to the field, then
% modify the selected region.  Matlab doesn't allow us access to the
% selection directly, so we'd have to use java. This is buggy and not
% forewards compatible, so disabled by default
useJavaAC = false;
try
	%	getting the java component that goes with the edit field so that we can
	%	set selection for doing tab completion later
	% this method of getting the java component probably doesn't work well in
	% newer versions of matlab; I should probably disable it
	j = get(get(h,'JavaFrame'),'FigurePanelContainer');
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
	drawnow;% to ensure the component exists before we try to find it
	f = j.getComponent(0).getComponent(0).getComponent(0);
	% if we did find the java component, set up to do autocomplete
	if(strcmp(f.getUIClassID,'TextFieldUI'))
		set(cs.name,'UserData',f);
	else
		useJavaAC = false; % this will be used by the callback function.  It's current value is captured when we create the anonymous functions below (I think; haven't really tested)
	end
catch
	useJavaAC = false;
end
set(cs.name,'KeyPressFcn',@nameAC,'Callback',@nameAC,'BusyAction','cancel');
% color selector
cs.color = uicontrol('background',[0 0 1]);
set(cs.color,'Callback',@(o,e)colorCB(o));
% marker type select
cs.marker = uicontrol( 'Style', 'popup',...
                       'String', ['x Cross|^ Upward-pointing triangle|Square|v Downward Pointing Triangle|'...
                            'o Circle|> Right-pointing triangle|Diamond|< Left-pointing triangle|'...
                            'Five-pointed star (pentagram)|Six-pointed star (hexagram)|. Point|+ Plus|* Asterisk'],...
                            ...%'*','^','s','v','o','>','d','<','p','h'}
                       'Callback', @(o,e)setMarkerCB(o));
marker = {'x','^','s','v','o','>','d','<','p','h','.','+','*'}; 
markerid = 0; markerset = []; markerList = get(cs.marker,'String');
% marker size selector
cs.markerSize = uicontrol( 'Style','edit','String','12','Callback',@(o,e)markerSizeCB(o),'backgroundcolor','w');
markerSize = 12;
% add button
cs.add = uicontrol('String','Add');
set(cs.add,'Callback',@(o,e)addCB(o))
% object list
cs.list = uicontrol('style','listbox','String',{'Feeder'},'UserData',true);
set(cs.list,'Callback',@(o,e)listCB(o));
% inspect button
cs.inspectB = uicontrol('String','Get Handles');
set(cs.inspectB,'Callback',@(o,e)disDelCB(o));
% enable/disable button
cs.enableB = uicontrol('String','Disable');
set(cs.enableB,'Callback',@(o,e)disDelCB(o));
% delete button
cs.deleteB = uicontrol('String','Delete');
set(cs.deleteB,'Callback',@(o,e)disDelCB(o))
% setup a callback for the data cursor
dc = datacursormode(h);
set(dc,'UpdateFcn',@(o,e)cursorCB(o,e));

% and update all their positions
set(h,'ResizeFcn',@(o,e)setControlPositions(o));
setControlPositions(cs);
drawnow;
circuit.cs = cs;
% finally update values for the the object list and name field
set(h,'UserData',circuit);
typeCB(cs.types);


%% perform data calculations

% build a struct lookup table for node coords
bl = circuit.buslist;
bl.id = fnSanitize(strtrim(cleanBus(bl.id)));
for i_=1:length(bl.id)
	bn = lower(['b' bl.id{i_}]);
	l.(bn).xy = bl.coord(i_,:);
	l.(bn).ind = i_;
	l.(bn).dev = {};
end

% build a hash lookup table for sections
sl = circuit.line;
buslist = lower(strcat('b',stripPhases({sl.bus1; sl.bus2}')));
mask = any(strcmp(buslist,'b'),2);
buslist = mat2cell(buslist,ones(1,size(buslist,1)),size(buslist,2));
buslist(mask) = [];
if(any(mask))
	namelist = sl(mask).Name;
	warning('circuitVisualizer:disconnectedSections','Some sections have no busses set! %s',sprintf('%s\n',namelist{:}));
end
namelist = strcat('line_',sl(~mask).Name);
% and, if applicatble, transformers
if(isfield(circuit,'transformer'))
	tl = circuit.transformer;
	namelist = vertcat(namelist,strcat('transformer_',{tl.Name}'));
	nbl = cell(length(tl),1);
	for i_=1:length(tl)
% 		nbl{i_} = tl(i_).Buses;
        nbl{i_} = stripPhases(tl(i_).Buses);
	end
	nbl = cellfun(@(x)strcat('b',x),nbl,'UniformOutput',false);
	buslist = [buslist; nbl];
end
namelist = fnSanitize(lower(namelist));
% actually add the lookup
for idx=1:length(buslist)
    buslist{idx} = strtrim(cleanBus(buslist{idx}));
end
buslist = cellfun(@fnSanitize,buslist,'UniformOutput',false);
for i_=1:length(buslist)
	bl = lower(buslist{i_});
	bn = namelist{i_};
	sec.(bn).bus = bl;
	sec.(bn).ind = i_;
	l.(bl{1}).dev{end+1} = bn;
	l.(bl{2}).dev{end+1} = bn;
end

% all the others
o_name = setdiff(circuit.n,{'buslist','bus','line','transformer'});
pmap = struct('load','bus1','generator','bus1','circuit','bus1','capacitor','Bus1');
missing_items = struct();
for objtype = o_name; objtype = objtype{1};
	ol = circuit.(objtype);
	if(~isfield(pmap,objtype)), continue; end;
	namelist = strcat([objtype,'_'],{ol.Name});
	bus_fn = pmap.(objtype);
	bl = strcat('b',stripPhases(ol.(bus_fn)));
	if(~iscell(bl)), bl = {bl}; end
	bl = lower(cleanBus(bl));
	for i_ = 1:length(ol)
		try % only mark buses that we know about
			l.(bl{i_}).dev{end+1} = namelist{i_};
		catch
			if(isempty('missing_items'))
				missing_items.name = namelist{i_};
			else
				missing_items(end+1).name = namelist{i_};
			end
		end
	end
end
if(~isempty(fieldnames(missing_items))), 
	warning('circuitVisualizer:missingItems','couldn''t locate bus to draw at for %i items',length(missing_items)); 
	for i_ = 1:length(missing_items)
		disp(missing_items(i_).name);
	end
end

% save the data for later use
circuit.bl = l;
circuit.sl = sec;

%% draw the whole circuit
% call the addToPlot subfunction to do the work
circuit.plots{1} = addToPlot(circuit, false, fieldnames(sec));

%% Finish creating
% make sure people don't overwrite the plot
set(h,'NextPlot','new');
% and save the data so we can get it back later by querying the figure
set(h,'UserData',circuit);

%% Object Selector Callback
	function typeCB(obj)
		% Here we handle new selections in the type and object selectors.
		% In particular, we make sure that when a new type is selected we
		% update the object list, and when either a new type or a new
		% object is selected, we update the name field
		
		c = get(get(obj,'Parent'),'UserData');
		cs = c.cs;
		% find the selected type
		type = c.n{get(cs.types,'Value')};
		% when type selector, we change the object selector as well
		if(obj==cs.types)
			if(strcmp(type,'bus')) %buslist data handled differently
				% notice that in either case we add an extra '[All]' item
				% to the object selector
				set(cs.objects,'String',[{'[All]'}; c.buslist.id]);
			else
				set(cs.objects,'String',[{'[All]'}; c.(type).Name]);
			end
			% default to select the "all" item
			set(cs.objects,'Value',1);
		end
		% get the selected object name
		oname = get(cs.objects,'String');
		if(iscell(oname))
			oname = oname{get(cs.objects,'Value')};
		end
		% set the name string
		% if we've selected "all" objects of a type, don't put that in the
		% name field
		if(strcmp(oname,'[All]'))
			set(cs.name,'String',type);
		else
			set(cs.name,'String',[type '.' oname]);
		end
	end

%% Name Autocompletion handler
	function nameAC(obj,e)
		% nameAC handles key events and callbacks from typing in the name
		% field.  This allows us to automatically select the correct items
		% in the type/object lists, and to insert the rest of the name in
		% the name field.
		% Unfortunately, doing autocomplete in the field requires accessing
		% the java objects directly in order to set the selected text, which is not consistent between matlab
		% versions and seems to lead to weird bugs anyway.
		
		% don't worry about special characters in keypressfcn calls
		if(isstruct(e) && ( isempty(e.Character) || (e.Character<' ' && ~strcmp(e.Key,'return')))), return; end
		
		% get the string and split it into type/object names
		try
			jo = get(obj,'UserData');
			oldt = char(jo.getText); % get(obj,'string') seems to be unreliable here for some reason...
		catch
			oldt = get(obj,'String');
		end
		n = regexp(oldt,'^(\w+)(\.)?(.*)?','once','tokens');
		if(isempty(n)), return; end; % can't work with strings that don't match
		
		% also collect the other objects that we'll need access to
		h = get(obj,'parent');
		c = get(h,'UserData');
		
		tl = get(c.cs.types,'String');
		if(isempty(n{2})) % no second part yet; autocomplete first part
			nind = strmatch(n{1},tl);
			if(isempty(nind)), return; end; % no matches
			newt = tl{nind(1)};
			newt(1:length(n{1})) = [];
			% go ahead and update the selected item in the type list and
			% the object list
			if(strcmp(tl{nind(1)},'bus'))
				ol = c.buslist.id;
			else
				ol = c.(tl{nind(1)}).Name;
			end
			if(nind(1)~=get(c.cs.types,'Value'))
				set(c.cs.types,'Value',nind(1));
				set(c.cs.objects,'String',[{'[All]'};ol],'Value',1);
			end
		elseif(~any(strcmp(n{1},tl)))
			% have a second part, but the first part doesn't match; can't
			% do anything
			return;
		else
			% we have a first part and it presumably matches, so make sure
			% we're showing the right list on the RHS, and get our list to
			% compare to
			nind = find(strcmp(n{1},tl),1);
			if(strcmp(n{1},'bus'))
				ol = c.buslist.id;
			else
				ol = c.(n{1}).Name;
			end
			if(nind~=get(c.cs.types,'Value'))
				set(c.cs.types,'Value',nind);
				if(get(c.cs.objects,'Value')>length(ol)+1)
					set(c.cs.object,'Value',1);
				end
				set(c.cs.objects,'String',[{'[All]'};ol]);
			end
			
			if(isempty(n{3}))
				% have a complete first part but can't work on second part
				% go ahead and pop up the list of values
				return;
			end
			% now start with autocomplete for the second part
			nind = strmatch(n{3},ol);
			if(isempty(nind)), return; end; % no matches
			if(iscell(ol))
				newt = ol{nind(1)};
			else
				newt = ol(nind(1),:);
			end
			newt(1:length(n{3})) = [];
			% if appropriate make sure the second column is properly
			% selected as well
			% with the java code enabled, it made some sense not to select
			% until we had a perfect match, but I think it works better to
			% just go for it.
			% nind = find(strcmp(n{3},ol),1);
			if(~isempty(nind))
				set(c.cs.objects,'Value',nind(1)+1);
			end
		end
		if(~isstruct(e) || isempty(newt)) % if we're doing a callback or don't have a string to add, don't autocomplete
			return;
		end;
		if(strcmp(e.Key,'return')) % go ahead and fill in the whole field permanently
			if(~isempty(n{3}))
				typeCB(c.cs.objects);
			else
				typeCB(c.cs.types);
			end
			return;
		end
		
		if(useJavaAC)
			% In java mode, add the new string to the text box and then
			% select the added text so that if the user keeps typing it
			% gets replaced
			
			% first set the string
			% nominally this causes the java object to update as well, but
			% in practice this seems to be unreliable.  Probably the
			% correct solution would be to make sure the subsequent java
			% calls are executed on the MATLAB Java EDT, but the version of
			% matlab we're working with doesn't support us using the EDT,
			% which is why I'm disabling java autocomplete by default.  If
			% you turn it on, this should work in a buggy way, or you can
			% uncomment the "set(jo,'text'..." line below, which works more
			% reliably but spews errors all over the console.
			set(obj,'String',[oldt newt]);
			%set(jo,'text',[oldt newt]);
			drawnow(); % force HG to update the java interface object so we can set selections
			jo.setSelectionStart(length(oldt));
			jo.setSelectionEnd(length(oldt)+length(newt));
		end
	end

%% Color selector callback
	function colorCB(obj)
		% pops up a color selector dialog and then sets the background
		% color of the requesting uicontrol.  Intended to allow the user to
		% easily select the color of the objects to be drawn
		color = uisetcolor(get(obj,'background'));
		set(obj,'background',color);
	end

%% Add button callback
	function addCB(obj)
		% This function is responsible for doing the UI stuff when adding
		% new data.  In particular, we determine which object(s) have been
		% requested, then get the data for them, and then use the addToPlot
		% function to actually plot the data.  Finally, we save the object
		% handles to apply changes (visibility, possibly style) to later
		%
		% The object data is determined by the currently selected item in
		% the type/object lists, but the name shown is from the name field.
		% Nominally these should match, but it's not a guarantee
		
		% grab data from the GUI objects
		c = get(get(obj,'Parent'),'UserData');
		cs = c.cs;
		oldlist = get(cs.list,'String');
		newlist = [oldlist;get(cs.name,'String')];
		% type name
		type = get(cs.types,'String');
		type = type{get(cs.types,'Value')};
		% object name
		name = get(cs.objects,'String');
		name = name{get(cs.objects,'Value')};
		
		% next we look up the data to add to the plot
		% Lines and transformers become line data, buses, loads, capacitors
		% and the like become point data, and some other objects (switches,
		% reclosers, etc) map to objects in one of the first two categories
		% in all cases, we will pass a line or bus name into the addToPlot
		% function, so we don't need to look up the actual data here.
		
		% first we resolve the references
		if(any(strcmp(type,{'switch','recloser','fuse'})))
			% switches, reclosers, and fuses refer to lines
			% lookup the line
			if(strcmp(name,'[All]'))
				name = c.(type).SwitchedObj;
			else
				if(~iscell(name)); name={name}; end;
				name = c.(type)(ismember(c.(type).Name,name)).SwitchedObj;
			end
			% remove the 'Line.' prefix...
			name = regexprep(name,'^Line\.','');
			% ... since we're recording that it's a line this way:
			type = 'line';
		elseif(strcmp(type,'linecode'))
			% linecodes map to a group of lines
			if(strcmp(name,'[All]'))
				name = c.(type).Name;
			elseif(~iscell(name)); name={name};
			end
			ind = ismember({c.line.LineCode},name);
			if(~sum(ind)), warndlg('No matching lines'), return; end
			name = c.line(ind).Name;
			type = 'line';
		end
		% now we handle the basic types
		switch(type)
			case {'line','transformer','regcontrol'}
				if(strcmp('regcontrol',type))
					type = 'transformer';
				end
				% these are line types
				if(strcmp(name,'[All]'))
					oname = c.(type).Name;
				else
					oname = name;
				end
				oname = strcat([type,'_'],oname);
				isNode = false;
			case 'bus'
				% buses are a special point type
				isNode = true;
				if(strcmp(name,'[All]'))
					oname = c.buslist.id;
				else
					oname = name;
				end
				oname = strcat('b',oname);
			case {'circuit','capacitor','capcontrol','storagecontroller','load','generator','vsource','pvsystem','storage','VviolationMax','VviolationMin'}
				% all these types are also point types that refer to a bus
				% however they refer to the bus by different names, so we
				% use this struct to determine what name to use to get the
				% bus
				pmap = struct('load','bus1','generator','bus1','circuit','bus1','vsource','bus1','capacitor','Bus1','capcontrol','Capacitor','storagecontroller','Element','pvsystem','bus1','storage','bus1','VviolationMax','Bus1','VviolationMin','Bus1');
				% capcontrol is special because it refers to a capacitor.
				% probably should've taken care of this in the previous
				% if/else block, but it works so I'm not gonna 'fix' it
				if(strcmp('capcontrol',type))
					if(strcmp(name,'[All]'))
						name = c.(type).(pmap.(type));
					else
						if(~iscell(name)); name={name}; end;
						name = c.(type)(ismember(c.(type).Name,name)).(pmap.(type));
					end
					type = 'capacitor';
				end
				if(strcmp('storagecontroller',type))
					if(strcmp(name,'[All]'))
						name = c.(type).(pmap.(type));
					else
						if(~iscell(name)); name={name}; end;
						name = c.(type)(ismember(c.(type).Name,name)).(pmap.(type));
					end
					type = 'storage';
				end
				isNode = true;
				if(strcmp(name,'[All]'))
					oname = c.(type).(pmap.(type));
				else
					if(~iscell(name)); name={name}; end;
					oname = c.(type)(ismember(c.(type).Name,name)).(pmap.(type));
				end
				oname = strcat('b',oname);
			otherwise
				% if you're getting this error, it might be as easy as
				% adding your new data type to one of the previous branches
				% of the switch statement.  Or you may need to write new
				% handling.  Or you may decide that you want to remove it
				% from the list, in which case you'll want to edit the code
				% where the type list is created (currently line 43 of this
				% file).
				errordlg('I don''t know how to graph that type of data yet!');
		end
		% once we've got an object name, we clean it up and add it to the
		% plot.  In particular, bus names may need .x.x.x suffixes that
		% specify phase connections to be removed.
		oname = stripPhases(oname);
		c.plots{length(newlist)} = addToPlot(c, isNode, oname, get(cs.color,'Background'),2);
		% update the GUI with the new plot list
		% the string value gets a new string, and the 'userdata' value
		% (which we're using to hold wheer each item is enabled or not)
		% gets a new 'true' at the end.
        newlist{end} = [num2str(size(oname,1)) ' ' newlist{end} ' - ' strtrim(markerList(markerid,:)) ' - ' get(cs.markerSize,'String') ' pts'];
		set(cs.list,'String',newlist,'Value',length(newlist));
        set(cs.list,'UserData',[get(cs.list,'UserData') true]);
		% and save the graph object handles so we can easily change their
		% properties later if we want
		set(get(obj,'Parent'),'UserData',c);
        legend('-DynamicLegend','Location','best');
	end

%% Plot List Callback
	function listCB(obj)
		% the main goal here is to set what text is shown on the
		% "enable/disable" button, depending on whether the selected item
		% is enabled or disabled.
		c = get(get(obj,'Parent'),'UserData');
		cs = c.cs;
		% recall that 'userdata' is a list with true/fales for whether the
		% corresponding list item is enabled or disabled.
		en = get(cs.list,'UserData');
		ind = get(cs.list,'value');
		en = en(ind);
		if(en)
			set(cs.enableB,'String','Disable');
		else
			set(cs.enableB,'String','Enable');
		end
	end

%% Callback for Delete and Dis/En-able buttons
	function disDelCB(obj)
		% delete, disable, or enable the selected plot items
		
		% first we get our data and figure out _which_ item we're talking
		% about
		c = get(get(obj,'Parent'),'UserData');
		cs = c.cs;
		i = get(cs.list,'Value');
		en = get(cs.list,'UserData');
		items = get(cs.list,'String');
		% Then we act on it
		if(obj==cs.enableB)
			% user clicked the enable/disable button.  Update the flag, the
			% button text, and the list item prefix, as well as setting or
			% unsetting the 'visible' property of the specified items.
			if(en(i))
				set(obj,'string','Enable');
				items{i} = ['(dis) ' items{i}];
				set(c.plots{i},'Visible','off');
			else
				set(obj,'String','Disable');
				items{i} = items{i}(7:end);
				set(c.plots{i},'Visible','on');
			end
			en(i) = ~en(i);
		elseif(obj==cs.inspectB)
			%inspect(c.plots{i});
			assignin('base','ans',c.plots{i});
		else % cs.deleteB
			% the user clicked the delete button.  Remove the item from the
			% list and delete the graph objects.
			en(i) = [];
			items(i) = [];
			if(i>length(en)); set(cs.list,'Value',length(en)); end;
			delete(c.plots{i});
			c.plots(i) = [];
			set(get(obj,'Parent'),'UserData',c)
		end
		% save the modified properties to the list data
		set(cs.list,'UserData',en);
		set(cs.list,'String',items);
	end

%% Control Positioning function
	function setControlPositions(cs)
		% this function is responsible for positioning the controls in the
		% window, which happens both at startup and any time the window is
		% resized.
		% We allow it to be called with either the struct of handles to the
		% controls or the handle of the figure.  We call the first way when
		% we setup the figure, and the resizefcn callback uses the second
		% method.
		if(ishandle(cs)) % called with the figure handle instead of the data we need
			cs = get(cs,'userdata');
			cs = cs.cs;
		end
		% get the window size for reference
		sz = get(get(cs.axes,'parent'),'position');
		sz(1:2) = [];
		% Set the positions of all the other controls
		% generally, the current design is a large graph at left.  At top
		% right, we have two selectors for type/object.  Underneath them is
		% a text box with the object's full name, along with a color picker
		% and a button to add it to the plot.  Underneath that is a list of
		% items currently in the plot, followed by buttons with which we
		% can hide or delete items from the list.
		set(cs.axes,'pos',[0 0 .65 1]);
		set(cs.types,'pos',[0.65*sz(1)+10 0.65*sz(2) sz(1)*0.15-15 0.35*sz(2)-10]);
		set(cs.objects,'pos',[0.8*sz(1)+5 0.65*sz(2) sz(1)*0.2-15 0.35*sz(2)-10]);
		set(cs.name,'pos',[0.65*sz(1)+10 0.65*sz(2)-30 sz(1)*0.15-15 25]);
        set(cs.color,'pos',[0.8*sz(1) 0.65*sz(2)-30+2 25-4 25-4]);
        set(cs.marker,'pos',[0.8*sz(1)+25 0.65*sz(2)-30 0.2*sz(1)-100 25-2]);
        set(cs.markerSize,'pos',[sz(1)-70 0.65*sz(2)-30+1 20 25-3]);
		set(cs.add,'pos',[sz(1)-47 0.65*sz(2)-30 40 25]);
		set(cs.list,'Position',[0.65*sz(1)+10 45 0.35*sz(1)-20 0.65*sz(2)-85]);
		set(cs.inspectB,'pos',[0.65*sz(1)+10 10 (0.35*sz(1)-40)/3 25]);
		set(cs.enableB,'pos',[0.7667*sz(1)+8 10 (0.35*sz(1)-40)/3 25]);
		set(cs.deleteB,'pos',[0.8833*sz(1)+5 10 (0.35*sz(1)-40)/3 25]);
	end

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
		oname = lower(fnSanitize(cleanBus(strtrim(oname))));
		
		% allocate an array to hold the handles for output.
		lo = zeros(length(oname),1);
		% lookup the data in the cached table and plot them one at a time
        counterrors=0;
        if isempty(markerset)
            markerset = 0;
        elseif ~markerset
            markerid = markerid + 1;
            set(cs.marker,'value',markerid);
            if markerid > length(marker), markerid = 1; end
        end
        newmarker = 1;
        if(isNode)
			for ii=1:length(oname)
				try
                    coords = c.bl.(oname{ii}).xy;
                    lo(ii) = plot(ah,coords(1),coords(2),marker{markerid},'Color',clr,'UserData',{isNode,oname{ii}},'LineWidth',line_width,'MarkerSize',markerSize);  
                    if newmarker
                        if length(oname) > 1
                            s = get(cs.types,'String');
                            s = s{get(cs.types,'value')};
                            s = [upper(s(1)) lower(s(2:end)) ' (' num2str(length(oname)) ')'];
                            set(lo(ii),'DisplayName',s);
                        else
                            if strcmpi(oname{ii},'bsourcebus')
                                set(lo(ii),'DisplayName', 'SourceBus (1)');
                            else
                                set(lo(ii),'DisplayName', [oname{ii} ' (1)']);
                            end
                        end
                        newmarker = 0;
                    else
                        set(lo(ii),'HandleVisibility','off');
                    end
                catch err
                    warning(err.message);
                    counterrors=counterrors+1;
                end
            end
            if counterrors>0; fprintf([num2str(counterrors) ' errors\n']);end
            else
			for ii=1:length(oname)
				try
					bn = c.sl.(oname{ii}).bus;
				catch
					fns = fieldnames(c.sl);
					[~, id] = ismember( lower(oname{ii}), lower(fieldnames(c.sl)) );
					bn = c.sl.(fns{id}).bus;
				end
				coords = vertcat(l.(bn{1}).xy,l.(bn{2}).xy);
				seg_l = hypot(diff(coords(:,1)),diff(coords(:,2)))/diag_l; % this should be the pixel length of this segment
				lo(ii) = plot(ah,coords(:,1),coords(:,2),'Color',clr,'UserData',{isNode,oname{ii}},'LineWidth',line_width,'MarkerSize',markerSize,'HandleVisibility','off');
				if(seg_l < 4)
					set(lo(ii),'Marker','x');
				end
			end
        end
        %
	end

%% Data Cursor Text
	function txt = cursorCB(o,e)
		% cursorCB is responsible for returning the custom datatip text
		% that shows what kind of object is selected.
		
		% t will be the target plot object
		t = get(e,'Target');
		% which in turn has a reference to it's data object and type
		oname = get(t,'UserData');
		isNode = oname{1};
		oname = oname{2};
		% we need the figure to get the circuit data from it
		h = get(get(t,'Parent'),'Parent');
		c = get(h,'UserData');
		% if something went very wrong and we can't find the circuit, just
		% use the object name
		if(~isstruct(c))
			txt = oname;
			return;
		end
		% otherwise, convert the name to a nice, readable form, then list
		% it along with whatever it's connected to
		% Note that each cell in our output becomes a line of the datatip
		if(isNode)
			txt = ['bus.' c.buslist.id{c.bl.(oname).ind}];
			txt = [{txt,'Devices:'}, regexprep(c.bl.(oname).dev,'([^_]+)_','$1.','once')];
		else
			txt = regexprep(oname,'([^_]+)_','$1.','once');
			bn = c.sl.(oname).bus;
			txt = [{txt,'Buses:'}, strcat('bus.',c.buslist.id(cellfun(@(x)c.bl.(x).ind,bn)))'];
			txt = [txt,'Devices:'];
			for i=1:length(bn)
				dl = c.bl.(bn{i}).dev;
				dl(strcmp(dl,oname)) = [];
				txt = [txt,regexprep(dl,'_','.','once')];
			end
		end
		
		% here we'll try to update the context menu for the datatip:
		try
			% get the datatip
			dt = e.DataTipHandle;
			if(isempty(dt)) % try to get it another way - the first works nicely in 2007; This seems to work in 2012a
				dt = get(datacursormode(h),'DataCursors');
				if(length(dt)>1)
					for i=1:length(dt)
						if(dt(i).Host==handle(e.Target))
							dt = dt(i); break;
						end
					end
				end
			end
			items = txt; items(~cellfun(@isempty,strfind(items,':'))) = [];
			% create a context menu
			m = uicontextmenu();
			uimenu(m,'Label','Add to Plot:','Enable','off');
			for i=1:length(items)
				uimenu(m,'Label',items{i},'Tag','1','Callback',@cursorMenuCB);
			end
			uimenu(m,'Label','Show Properties:','Enable','off');
			for i=1:length(items)
				uimenu(m,'Label',items{i},'Tag','0','Callback',@cursorMenuCB);
			end
			dt.UiContextMenu = m;
		catch
			disp(75);
		end
	end

%% Callback for cursor context menu
	function cursorMenuCB(o,e)
		n = regexp(get(o,'Label'),'^(\w+)\.(.*)','once','tokens');
		c = get(get(get(o,'Parent'),'Parent'),'UserData');
		if(strcmp(get(o,'Tag'),'1')) % add to plot
			tind = find(strcmp(n{1},get(c.cs.types,'String')));
			set(c.cs.types,'Value',tind);
			typeCB(c.cs.types); % update the object list
			oind = find(strcmp(n{2},get(c.cs.objects,'String')));
			set(c.cs.objects,'Value',oind);
			typeCB(c.cs.objects);
			addCB(c.cs.add);
		else
			if(strcmp(n{1},'bus'))
				b = c.bl.(['b' fnSanitize(n{2})]);
				fprintf('Bus %-13s (%d,%d)\n',[n{2} ':'],b.xy(1),b.xy(2));
				disp(regexprep(b.dev(:),'_','.','once'));
			else
				oind = strcmp(c.(n{1}).Name,n{2});
				display(c.(n{1})(oind));
			end
		end
    end

%% callback for setting marker
    function setMarkerCB(o,e)
        markerset = 1;
        markerid = get(o,'Value');
    end

%% callback for setting marker
    function markerSizeCB(o,e)
        markerSize = str2double(get(o,'String'));
    end
end
