function a = subsasgn(a,index,val)
% SUBSASGN Define index assignment for asset objects

% Sbusasgn is a tricky function because of the way matlab handles most
% data: instead of working with pointers to objects, when you assign an
% object value to a new variable, it makes a lazy copy of the object, which
% is detached from the original as soon as you assign anything to it.  So
% the strategy that we might normally do of "follow the reference down the
% chain and then assign at the last level" isn't good enough, because it
% would modify a copied object and then throw away the changes when the
% function returns.
%
% Instead, we work down the chain to look up the right object to assign to,
% and then work back up the chain assigning in the modified results.
% Fortunately, unless there are bugs, you shouldn't need to modify this.
%
% newer versions of matlab are nice in that the classdef syntax takes care
% of a lot of this for you (never tried overriding subsasgn as well), and
% handle classes work a little bit more like pointers.  But we're stuck in
% R2007b, so whatever...

%% Work down the tree to find a value to assign to
% fill in "values we're going to modify" in the v cell array
v{1} = a;
v{length(index)+1} = val;
for i = 2:(length(index))
	switch(index(i-1).type)
		case '()' % array type indexing
			% be careful when we might be dynamically resizing the array!
			try
				v{i} = v{i-1}(index(i-1).subs{:});
			catch err
				if(~strcmp(err.identifier,'MATLAB:badsubscript')), rethrow(err); end
				v{i-1}(index(i-1).subs{:}) = eval(mfilename('class'));
			end
		case '.' % struct style indexing
			% If it's our class, get the subproperty using the "get"
			% function.  If not, rely on whatever else it is to handle the
			% subsasgn for us and short-circuit the rest of the steps
			if(isa(v{i-1},mfilename('class')))
				v{i} = get(v{i-1},index(i-1).subs);
			else
				v{i-1} = subsasgn(v{i-1},index(i-1:end),val);
				index(i-1:end) = [];
				break;
			end
	end
end
%% Work back up, assigning the modified values
% really makes you miss pointers...
for i=length(index):-1:1
	switch index(i).type
		case '()' % array style indexing
			%if we're trying to turn an empty variable into a class, matlab
			%will complain, so we explicitly assign it to be of the class
			%first before assigning to the desired index
			if(isempty(v{i}))
				v{i} = eval(mfilename('class'));
			end
			if(isempty(v{i+1})) %using assignment to reduce size
				v{i}(index(i).subs{:}) = [];
			else
				v{i}(index(i).subs{:}) = v{i+1};
			end
		case '.' % struct style indexing
			if(isa(v{i},mfilename('class')))
				% for members of our class, only set if we recognize the
				% field name, and if so, then just call set() to do the
				% hard work
				fn = fieldnamefix(index(i).subs,a(1).fieldnames);
				if(~isempty(fn))
					v{i} = set(v{i},fn,v{i+1});
				else
					error('unknown property, ''%s''',index(i).subs);
				end
			else
				% once in a while we get a struct or some other class that
				% is the object that we're assigning the last value into,
				% so it hasn't broken out in the 'else' case on the way
				% down the loop; again we just cal subsasgn on that object
				% and let it figure out what to do
				v{i} = subsasgn(v{i},index(i),v{i+1});
			end
	end
end
% don't forget to write the results back to the input variable when we're
% all done so that our changes are saved!
a=v{1};

end
