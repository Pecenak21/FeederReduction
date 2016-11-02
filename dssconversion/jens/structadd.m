function ns = structadd(s_source,fn_source,fn_source_ID,s_dest,fn_dest,fn_dest_ID,NewID,fn_NewID)
% structadd(s_source,fn_source,fn_source_ID,s_dest,fn_dest,fn_dest_ID,NewID,fn_NewID)
%
% PURPOSE : Add information contained in s_source.fn_source to
% s_dest.fn_dest. All fields that already exist in fn_dest and that
% do not exist in fn_source will be padded with ''. All fields that exist
% in fn_source, but not in fn_dest will be added to fn_dest and padded with
% ''.
%
%
% INPUT :   s_source: struct that has field with info to be added
%           fn_source: fieldname of field with info
%           fn_source_ID: ID of source
%           s_dest: struct that has field where new info will be added
%           fn_dest: fieldname of field where new info will be added
%           NewID: Identifier that specifies type of new object
%           fn_NewID: fieldname where NewID is specified
%
% OUTPUT :  a structure that contains information from the destination
% structure and from the source structure. Uniform number of entries is
% achieved by padding with ''.

if ~isfield(s_dest,fn_dest)
    s_dest.(fn_dest).(fn_source_ID)=''; % destination does not exist, yet => initialize
end
if isfield(s_source,fn_source) 
    fn=GetFieldNames('unique',s_dest,fn_dest,s_source,fn_source);
    n_source=length(s_source.(fn_source).(fn_source_ID));
    if isfield(s_dest.(fn_dest),fn_dest_ID)
        n_dest=length(s_dest.(fn_dest).(fn_dest_ID));
    else
        n_dest=0;
    end
    for i=1:length(fn)
        if isfield(s_dest.(fn_dest),fn{i}) && isfield(s_source.(fn_source),fn{i}) % field exists in both structs, just add content in .switch to .structure
                s_dest.(fn_dest).(fn{i})=[s_dest.(fn_dest).(fn{i});s_source.(fn_source).(fn{i})]; 
        elseif isfield(s_source.(fn_source),fn{i}) % field only exist in .switch, need to pad beginning of .(fn_dest)
            clear dummy;
            
            if ischar(s_source.(fn_source).(fn{i}){1})
                dummy(1:n_dest,1)=deal({''});
            else
                dummy(1:n_dest,1)=deal({[]});
            end
            s_dest.(fn_dest).(fn{i})=[dummy;s_source.(fn_source).(fn{i})];
        elseif isfield(s_dest.(fn_dest),fn{i}) % field only exist in .(fn_dest), need to pad end of .(fn_dest)
            if strcmp(fn{i},fn_NewID)
                clear (fn_NewID);
                Type(1:n_source,1)=deal({NewID});
                s_dest.(fn_dest).(fn{i})=[s_dest.(fn_dest).(fn{i});Type];
            else
                clear dummy;
                if ischar(s_dest.(fn_dest).(fn{i}){1})
                    dummy(1:n_source,1)=deal({''});
                else
                    dummy(1:n_source,1)=deal({[]});
                end
%                 dummy(1:n_source,1)=deal({''});
                s_dest.(fn_dest).(fn{i})=[s_dest.(fn_dest).(fn{i});dummy];
            end
        end
        if ~isfield(s_dest.(fn_dest),fn_NewID)
            clear (fn_NewID);
            Type(1:n_source,1)=deal({NewID});
            s_dest.(fn_dest).(fn_NewID)=Type;
        end
    end
end
ns=s_dest;

end
