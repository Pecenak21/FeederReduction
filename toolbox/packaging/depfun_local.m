function list = depfun_local(file, basepath, looked_up)
% depfun_local(file,[basepath]) returns a list of dependent functions, only including files in a specific directory.
% it does not currently attempt to duplicate any of the other outputs of depfun.
%
% if no basepath is specified, the directory in which _file_ resides will be used.
%
% the third argument is for internal use only

% implement the default basepath
if(nargin<2)
	basepath = fileparts(which(file));
end
if(nargin<3)
	looked_up = {};
end

% get top-level dependencies
looked_up{end+1,1} = which(file);
list = depfun(file,'-toponly','-quiet');
% remove things outside the base path
m = cellfun(@isempty,regexp(list,['^' basepath],'once'));
list(m) = [];

% recurse on paths that we still need to look up
todo = list(~ismember(list,looked_up));
while(~isempty(todo))
	% get the list of new dependent functions, and add them to our master list, as well as the list of functions we don't need to check anymore
	newlist = depfun_local(todo{1},basepath,looked_up);
	looked_up = vertcat(newlist,looked_up);
	list = vertcat(list,newlist);
	
	% refresh the todo list
	todo = list(~ismember(list,looked_up));
end

% probably have some duplicates, since the recursive calls return the function they were called on, which is already on the list
list = unique(list);

end
