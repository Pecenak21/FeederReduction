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
		case {'name','like','bus1','bus2'}
            val{i} = dataclean(val{i},'name');
		case {'linecode','geometry'}
            if ischar(val{i}) % if name of the object is passed in
                val{i} = dataclean(val{i},'name');
            else % if an dss object is passed in
                % get fieldnames for argument object
                fns1 = fieldnames(val{i});
                % get fieldnames for target object
                fns2 = fieldnames(s);
                % search for matching fields between target and argument objects
                [id1 id2] = ismember( lower(fns1) , lower(fns2) );
                for ii = 1:length(id1)
                    if id1(ii) && ~isempty( val{i}.(fns1{ii}) )
                        s.defaults.(fns2{id2(ii)}) = val{i}.(fns1{ii});
                    end
                end
                val{i} = val{i}.Name;
            end
            if ~isempty(s.data.LineCode) && ~isempty(s.data.Geometry)
                warning('dssline:multiplespecification','You have specified both LineCode and LineGeometry, the later will overwrite the previous one. Be careful!!!');
            end
        case 'phases'
			val{i} = dataclean(val{i},'phase'); 
		% Emergamps defaults to 1.35x normamps
		case 'normamps'
			val{i} = dataclean(val{i},'num');
			s.defaults.Emergamps = val{i}*1.35;
		case {'length','basefreq','emergamps','faultrate','pctperm','rg','xg','rho','repair','c1','c0','cmatrix'}
			val{i} = dataclean(val{i},'num');
		case {'r1','x1','r0','x0','rmatrix','xmatrix'}
			val{i} = dataclean(val{i},'num');
% 			if(get(s,'Phases')==3)
				s.data.(propn{i}) = val{i}; % preemptively set so we can use get() below
				if(length(propn{i})==2) % sequence form; calc matrix values
					z = get(s,{'R0','X0','R1','X1'});
					[s.defaults.Rmatrix, s.defaults.Xmatrix] = zconv(z{:});
					if(~isempty(s.data.Rmatrix)||~isempty(s.data.Xmatrix))
						warning('dssline:impedanceoverride','Overriding specified impedance matrix with sequence impedances');
						s.data.Rmatrix = [];
						s.data.Xmatrix = [];
					end
				else % matrix form; calc sequence values
					z = get(s,{'Rmatrix','Xmatrix'});
					% go ahead and clear the data values
					s.data.R0 = []; s.data.X0 = [];
					s.data.R1 = []; s.data.X1 = [];
					me = warning('error','zconv:sketchyresult');
					try
						[s.defaults.R0, s.defaults.X0, s.defaults.R1, s.defaults.X1] = zconv(z{:});
					catch
						[s.defaults.R0, s.defaults.X0, s.defaults.R1, s.defaults.X1] = deal(nan);
					end
					warning(me);
				end
% 			else
% 				disp(char(s));
% 				warning('I haven''t learned yet how to use sequence impedances for lines with phases != 3');
% 			end
		case {'switch','enabled'}
			val{i} = dataclean(val{i},'logical','string');
        case 'units'
			sunits = get(s,'Units');
			% if the new units don't match the old units, we need to convert
            if(~(strcmpi(sunits,'none') || strcmpi(sunits,val{i})))
				% For converting units, our strategy is to come up with a 'uconv' parameter 
                % which allows us to convert the old length via multiplication (i.e. length = length*uconv;), 
                % which we will also use to convert all the impedance quantities
				% In generating uconv, we first go to meters, and then on to our destination unit
				uconv = 1;
				switch(lower(sunits(1)))
					case 'm'
						uconv = 1;
					case 'f'
						uconv = 12*2.54/100;
					case 'k'
						if(lower(sunits(2)) == 'f') % kft
							uconv = 12*2.54/100*1000;
						else %km
							uconv = 1000;
						end
				end
				switch(lower(val{i}(1)))
					case 'm'
						% uconv = uconv;
					case 'f'
						uconv = uconv/(12*2.54/100);
					case 'k'
						if(lower(sunits(2)) == 'f') % kft
							uconv = uconv/(12*2.54/100*1000);
						else %km
							uconv = uconv/1000;
						end
					otherwise
						warning('dssline:units','trying to convert to units but can''t understand your new units. they will be interpreted as meters.');
				end
				% convert lengths
				s.data.Length = s.data.Length*uconv; s.defaults.Length = s.defaults.Length*uconv;
				% convert impedances
				for fn = {'R1','R0','X1','X0','C1','C0','Rmatrix','Xmatrix','Cmatrix','Rg','Xg'}
					s.data.(fn{1}) = s.data.(fn{1})/uconv;
					s.defaults.(fn{1}) = s.defaults.(fn{1})/uconv;
				end
			end
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
	s.data.(propn{i}) = val{i};
end

end
