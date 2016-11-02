function ns = structreduce(struc,index,flgUnique)
% structreduce only keeps the objects in a structure that are included in
% the index array, index can be made unique by setting flgUnique

if nargin<3
    flgUnique=1;
end

% get the field names
fn = fieldnames(struc)';

if flgUnique % only unique indeces
    index=unique(index);
    
   
end
% no zero indeces
    index=index(index~=0);
% index=find(index);

if(length(struc)==1)
	% extract data to go with the fieldnames
	for i=1:length(fn);
        struc.(fn{i})= struc.(fn{i})(index);
	end
else
	error('Only works with structures that have size "one"')
end

% pass the fieldnames and data cells into struct() to get a new struct
ns = struc;

end
