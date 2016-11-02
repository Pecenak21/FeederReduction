function saveFile(filepath,dat,overwrite)
global indent; if isempty(indent), indent = ''; end
if ~exist('filepath','file') || (exist('overwrite','var') && ~overwrite)
    try save(filepath,'-struct','dat');
    catch e
        save(filepath,'dat');
    end
    fprintf(['%sSaved file: ' filepath '\n'],indent);
elseif exist('filepath','file')
    fprintf(['%sFile exists. Skip! File path: ' filepath '\n'],indent);
end
end