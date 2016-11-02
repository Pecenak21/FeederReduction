function ns = SectionTrim(s,fn_ID,ID_delete,fn_b1,fn_b2,fn_b3,fn_closed)
% SectionTrim(s,fn_ID,ID,fn_b2,fn_b3)
%
% PURPOSE : Delete a section and replace the bus id that connects the
% deleted section with another section with the bus id of the other side of
% the deleted section.
%
%
% INPUT :   s: struct that has all entries (e.g., o.Section)
%           fn_ID: fn that identifies deleted section (e.g.'SectionType')
%           ID_delete: ID of section that will be deleted (e.g. 'switch')
%           fn_b2: field name with common bus ID (e.g., ToNodeID) 
%           fn_b3: field name with bus ID on the other end of section
%
% OUTPUT :  
% before:
% b1---line/cap/transformer---b2--switch--b3--line/cap/transformer--b4
% or
% b2--line/cap/transformer--b1--switch--b3--line/cap/transformer--b4
%
% after: 
% b1--line/cap/transformer--b3--line/cap/transformer--b4

if nargin<7
    flgAllClosed=1;
else
    flgAllClosed=0;
end

% get indeces of object that will be deleted
flg=ismember(s.(fn_ID),ID_delete);
index=find(flg);

% if SectionID of object to be deleted is a SectionID of object to be
% retained, the simply delete the former (this means switch and line are in
% parallel)
SectionID_switch=s.SectionID(flg);
fn=fieldnames(s);
n_fn=length(fn);
for i=1:length(SectionID_switch)
    if ~ismember(SectionID_switch(i),s.SectionID(~flg))
        flg(index(i))=0;
    end
end
% remove fields associated with ID_delete
fn=fieldnames(s);
n_fn=length(fn);
for i=1:n_fn
    s.(fn{i})(flg)=[];
end


% Deal with the remaining switches, the code below is not verified
flg=ismember(s.(fn_ID),ID_delete);
if any(flg)
    % b2 in s will be replaced by b3
    b2=s.(fn_b2)(flg);
    b3=s.(fn_b3)(flg);
    if ~flgAllClosed
        section_closed=s.(fn_closed)(flg);
    end


    % remove fields associated with ID_delete
    fn=fieldnames(s);
    n_fn=length(fn);
    for i=1:n_fn
        s.(fn{i})(flg)=[];
    end

    % replace b2 with b3
    b=[b2,b3];
    n_b=length(b2);
    for i=1:n_b
        flg=strcmp(s.(fn_b2),b2{i});
        if any(flg) && ~flgAllClosed
            if strcmp(section_closed{i},'NONE')
                s.(fn_b2)(flg)=b(i,2);
            end
        else
            s.(fn_b1)(flg)=b(i,2);
        end
        flg=strcmp(s.(fn_b1),b2{i});
        if any(flg) && ~flgAllClosed
            if strcmp(section_closed{i},'NONE')
                s.(fn_b1)(flg)=b(i,2);
            end
        else
            s.(fn_b1)(flg)=b(i,2);
        end
    end
end

ns=s;

end
