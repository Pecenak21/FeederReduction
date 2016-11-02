function c = rmBusByName(c,bname)
id = ismember(lower(c.buslist.id),lower(bname));
if ~isempty(id)
    c.buslist.id(id) = [];
    c.buslist.coord(id,:) = [];
end
end