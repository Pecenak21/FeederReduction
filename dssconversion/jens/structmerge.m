function ns = structmerge(struct1,ID1,struct2,ID2,flgDeleteEmpty,flgStruct2)
% Structmerge finds the entries with commone IDs in struct1 and struct2
% and adds the corresponding objects in struct2 to struct1.
% If flgDeleteEmpty is set, then fields that only have '' entries are
% removed.

if nargin<5
    flgDeleteEmpty=0;
    flgStruct2=0;
end

[flg,index]=ismember(struct1.(ID1),struct2.(ID2));

if all(index)
    fn=fieldnames(struct2);
    for i=1:length(fn)
        if ~isfield(struct1,fn{i}) % only add field if it does not exist already
            struct1.(fn{i})=struct2.(fn{i})(index);
        elseif flgStruct2
            struct1.(fn{i})=struct2.(fn{i})(index); % struct2 replaces stuff in struct1
        end
    end
end

if flgDeleteEmpty
    fn=fieldnames(struct1);
    for i=1:length(fn)
        if ~isnumeric(struct1.(fn{i}){1}) % keep field if it is numeric
            [flg,index]=ismember( struct1.(fn{i}),'');
            if all(index)
                struct1=rmfield(struct1,fn{i});
            end
        end
    end
end


ns=struct1;
end
