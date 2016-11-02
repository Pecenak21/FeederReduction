function [d] = getGmatrix(c,comp,bus0)
global indent; if isempty(indent), indent = ''; end

if strfind(class(comp),'dss')
    g = comp; type = strrep(class(comp),'dss','');
elseif ischar(comp)
    g = c.(comp); type = lower(g);
else
    error('not supported ''comp'' input!');
end

fp = fNamePrefix(['Gmat_' type]); fp = [fp{1} '.mat'];
if exist(fp,'file')
    fprintf(['%sGmatrix saved file exists: ' fp '. Load to use.\n'],indent);
    d = load(fp); d = d.gmat;
    return;
end
comp2disable = {'generator','pvsystem','storage'};
id = find(ismember(comp2disable,fieldnames(c))); 
if ~isempty(id)
    for i = 1:id
        for k = 1:length(c.(comp2disable{id(i)})), c.(comp2disable{id(i)})(k).enabled = 'no'; end
    end
end
d.gmat = calculateGMatrix([{g.bus1} bus0],c);
saveFile(fp,d);
fprintf(['%sSaved Gmatrix file: ' fp '\n'],indent);
end