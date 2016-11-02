function ns = structconv(os)
% structconv converts from a struct array to a struct whose values are cell
% arrays or vice versa.

% get the field names
if ~isempty(os)
    cellform = fieldnames(os)';


    if(length(os)==1) 
        if ischar(os.(cellform{1}))
            data = struct2cell(os);
            % reshape the data into the new cell array
            for i=1:length(cellform);
                cellform{2,i} = {data(i,:)'};
            end
        else
            % extract data to go with the fieldnames
            for i=1:length(cellform);
                cellform{2,i} = os.(cellform{1,i});
            end
        end
    else
        data = struct2cell(os);
        % reshape the data into the new cell array
        for i=1:length(cellform);
            cellform{2,i} = {data(i,:)'};
        end
    end

    % pass the fieldnames and data cells into struct() to get a new struct
    ns = struct(cellform{:});
else
    ns=[];
end

    
