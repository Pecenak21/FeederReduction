function saveFigure(fhandle,filepath,overwrite)
global indent; if isempty(indent), indent = ''; end
if ~exist('filepath','file') || (exist('overwrite','var') && ~overwrite)
    saveas(fhandle,filepath);
    disp([indent 'Saved file: ' filepath ]);
elseif exist('filepath','file')
    disp([indent 'File exists. Skip! File path: ' filepath ]);
end
end