function s = set(s,propn,val,varargin)
% The Set function is responsible for handling all the magic of making sure
% we format the object data correctly.  For example:
%	* Some kinds of data we try to correct subtle mistakes
%	* Setting some properties in opendss triggers other properties to be
%	cleared or have a default value
% All of this kind of detail is taken care of here in the set function

% handle cell value
while iscell(val) && length(val)==1
	val = val{1};
end
while iscell(propn) && length(propn)==1
	propn = propn{1};
end

% handle setting multiple values at once
if(nargin > 3)
	val = reshape([{propn,val} varargin],[2 (nargin-1)/2])';
	propn = fieldnamefix(val(:,1),s.fieldnames);
	val = val(:,2);
else %otherwise just make sure we're setting a valid value
	propn = fieldnamefix({propn},s.fieldnames);
	val = {val};
end

% set the data, altering a few cases as desired and setting other values or
% defaults in some other cases
for i=1:length(propn)
	if(isempty(propn{i})), continue; end
	% as we go through the main body of the loop, we'll tweak the value
	% itself, and make any necessary changes to OTHER properties, and then
	% let one line be responsible for setting the current property
	switch(lower(propn{i}))
		% some programs specify phases individually, but we just want the number
		case 'phases'
			val{i} = dataclean(val{i},'phase');
		% EmerghkVA defaults to 1.5x normhkVA
		case 'normhkva'
			val{i} = dataclean(val{i},'num');
			s.defaults.EmergHKVA = val{i}*1.5;
		case {'kv','r','rneut','xneut','maxtap','mintap','numtaps'}
			val{i} = dataclean(val{i},'num');
		case 'tap'
			val{i} = dataclean(val{i},'num');
		case 'kva'
			val{i} = dataclean(val{i},'num');
			if(get(s,'Wdg')==1)
				% when assigning a kva rating for winding 1, OpenDSS makes
				% that the value for all other windings.  Here, we take a
				% weaker form and don't overwrite already assigned values
				s.defaults.kVA(:) = val(i);
				s.defaults.kVAs(:) = val{i};
			end
        case 'windings'
			val{i} = dataclean(val{i},'num');
			windingexpand(val{i},false);
            s.defaults.Wdg = 1;
		case 'wdg'
			val{i} = dataclean(val{i},'num');
			% expand the total number of windings if necessary
			if(val{i}>get(s,'Windings'))
				windingexpand(val{i});
			end
			if(val{i}<1), val{i} = 1; end
        case {'conns','conn'}
			[val{i} grounded] = dataclean(val{i},'conn');
			if(length(grounded) > get(s,'Windings'))
				windingexpand(length(grounded));
			end
			
			% if a terminal is grounded, update bus from default ungrounded
			% bus to grounded bus (.1.2.3.0)
			if(iscell(grounded)), grounded = [grounded{:}]; end
			for j = find(grounded)
% 				if( get(s,'Phases')==3 && ~isempty(s.data.Buses{j}) )
% 					s.data.Buses{j} = [s.data.Buses{j} '.1.2.3.0']; % js: commented out for now as it results in something like that .1.2.3.1.2.3.0
% 				elseif(strcmp(s.data.Buses{j}(end-1:end),'.0')) % implied phases = 1
% 					warning('dsstransformer:grounding','''%s'' may need to be explicitly grounded',s.data.Name);
% 				end
			end
		case {'buses','bus'}
			val{i} = dataclean(val{i},'name');
			if(iscell(val{i}) && length(val{i})>get(s,'Windings'))
				windingexpand(length(val{i}));
			end
        case 'name'
			val{i} = dataclean(val{i},'name');
			% update default bus names
			s.defaults.Bus =  strcat(val{i}, '_',cellfun(@num2str,num2cell(1:get(s,'Windings')),'UniformOutput',false));
			s.defaults.Buses = s.defaults.Bus;
		case 'like'
			if(ischar(val{i}))
				val{i} = dataclean(val{i},'name');
			elseif(isa(val{i},mfilename('class')))
				s.defaults = get(val{i});
				val{i} = s.defaults.Name;
            end
        case {'enabled'}
			val{i} = dataclean(val{i},'logical','string');
		otherwise
			% determine property's type/class and convert given value to that type
			cl = class(s.defaults.(propn{i}));
			if ~isa(val{i},cl)
				switch cl
					case 'char'
						val{i} = char(val{i});
					case 'double'
						val{i} = dataclean(val{i},'num');
				end
			end
	end
	pla = {'Buses','Conns','kVs','kVAs','Taps','Rs'};
	plb = {'Bus','Conn','kV','kVA','tap','R'};
	% Special treatment for multiple-valued properties
	try s.dependency; hasdep=true; catch, hasdep=false; end
	if(hasdep && ismember(propn{i},s.dependency(2:end)))
		% lookup the index we're currently using
		id = get(s,s.dependency{1});
		% assign the data
		s.data.(propn{i}){id} = val{i};
		% make the values show up on the plural properties
		if(~ismember(propn{i},plb)), continue; end
		fidx = find(strcmp(propn{i},plb));
		% need to remove the old plural data so the new stuff can show
		% through, but first we have to save it in:
		%	* plural defaults;  * multi-value data
		if(~isempty(s.data.(pla{fidx})))
			olddata = s.data.(pla{fidx});
			s.data.(pla{fidx}) = [];
			s.defaults.(pla{fidx})(1:length(olddata)) = olddata;
			if(fidx>2), olddata = num2cell(olddata); end
			s.data.(propn{i})(1:length(olddata)) = olddata;
			% reassign the new data 'cause I can't think of a better way
			s.data.(propn{i}){id} = val{i};
		end
		if(fidx<=2) %Bus/Con - i.e. the string types
			s.defaults.(pla{fidx}){id} = val{i};
		else % numeric types
			s.defaults.(pla{fidx})(id) = val{i};
		end
	else
		s.data.(propn{i}) = val{i};
	end
	% And for the plural forms of properties, apply them to the
	% multiple-valued properties
	[ism_ ism_] = ismember(propn{i},pla);
	if(ism_)
		if(~iscell(val{i}))
			if(ischar(val{i}))
				val{i} = val(i); %wrap char strings as cells
			else
				val{i} = num2cell(val{i}); % everything else gets converted from an array to a cell array
			end
		end
		% set our data as defaults, so it'll show up but not be written to
		% file
		s.defaults.(plb{ism_})(1:length(val{i})) = val{i};
		% blank the data elements so the defaults can show through
		s.data.(plb{ism_})(1:length(val{i})) = {[]};
	end
end

	function windingexpand(num,dowarn)
		if(~exist('dowarn','var') || isempty(dowarn) || dowarn)
			warning('dsstransformer:addWinding','Automatically adding one or more windings to ''%s''.  This is not how opendss behaves, so if you got this warning parsing a dss file, you may get different results',s.data.Name);
		end
		s.data.Windings = num; % we extend 'data' not defaults because opendss doesn't do this automatically, so we _need_ to write that to file
		s.defaults.Bus(end+1:num) =  strcat(get(s,'Name'),cellfun(@num2str,num2cell(length(s.defaults.Bus)+1:num),'UniformOutput',false));
		% expand the data and defaults arrays for the all the multi-valued
		% properties
		for j = 2:length(s.dependency)
			s.data.(s.dependency{j}){num} = [];
			s.defaults.(s.dependency{j})(end+1:num) = s.defaults.(s.dependency{j})(end);
		end
		% expand only the defaults arrays for the plural properties
		oldl = length(s.defaults.Buses);
		s.defaults.Buses(oldl+1:num) = s.defaults.Bus(oldl+1:num);
		s.defaults.Conns(oldl+1:num) = s.defaults.Conn(oldl+1:num);

		fns_ = {'kVs','kVAs','Taps','Rs';...
				'kV', 'kVA', 'tap', 'R'};
		for j=i:length(fns_);
			s.defaults.(fns_{1,j})(oldl+1:num) = [s.defaults.(fns_{2,j}){oldl+1:num}];
		end
	end

end
