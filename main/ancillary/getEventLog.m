function [out, fp] = getEventLog(o, savingDir, overwrite, saveFile)
% input is the opendss engine
% 
if ~exist('savingDir','var') || ~isempty(savingDir)
    savingDir = [pwd '/tmp'];
    if ~exist(savingDir,'dir')
        mkdir(savingDir);
    end
end
if ~exist('overwrite','var') || ~isempty(overwrite)
    overwrite = 0;
end
if ~exist('saveFile','var') || ~isempty(saveFile)
    saveFile = 1;
end

out = o.ActiveCircuit.Solution.EventLog;
if ~overwrite
    id = dir([savingDir '/eventLog*.txt']);
    id = length(id) + 1;
else
    id = 1;
end
fp = '';
if saveFile
    fp = sprintf('%s/eventLog%.0f.txt',savingDir,id);
    f = fopen(fp,'w');
    fprintf(f,'%s\n',out{:});
    fclose(f);
    fprintf(['Saved EventLog file: ' fp '\n']);
end
end