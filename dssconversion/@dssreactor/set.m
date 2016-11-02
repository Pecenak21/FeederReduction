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
		% Emergamps defaults to 1.35x normamps
		case 'normamps'
			val{i} = dataclean(val{i},'num');
			s.defaults.Emergamps = val{i}*1.35;
		% standardize connection type names
		case {'kvar','kv','basefreq','r','rmatrix','rp','x','xmatrix','emergamps','faultrate','pctperm','repair'}
			val{i} = dataclean(val{i},'num');
		case {'parallel','enabled'}
			val{i} = dataclean(val{i},'logical','string');
		case 'conn'
			[val{i} grounded] = dataclean(val{i},'conn');
			
			% In case the capacitor is not grounded, we need to specify the
			% second bus as .4.4.4.
			if ~grounded
				s.data.Bus2 = [s.data.Bus1 '.4.4.4'];
			end
		case {'name','like','bus1','bus2'}
			val{i} = dataclean(val{i},'name');
		otherwise
			% determine property's type/class and convert given value to that type
			warning('dssreactor:unknownparameter',['parameter' propn{i} ' doesn''t exist. check again. u shouldn''t end up here']);
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
