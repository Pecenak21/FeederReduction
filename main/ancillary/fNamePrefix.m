function fn = fNamePrefix(prefixStr,r,outDir)
% generate prefix name with output folder path in it for saving simulation's results or figures
if ~exist('prefixStr','var') || isempty(prefixStr)
    prefixStr = '';
end
if ~exist('outDir','var') || isempty(outDir)
    outDir = [r.conf.outputDir '/' prefixStr];
else
    outDir = [outDir '/' prefixStr];
end
if length(r) == 1
    fn = {fNamePrefixElmt(prefixStr,r,outDir)};
else
    for i = 1:length(r)
        fn{i,1} = fNamePrefixElmt(prefixStr,r(i),outDir);
    end
end
end
function fn = fNamePrefixElmt(prefixStr,r,outDir)
global fName; global fdSetup; global tDay; global pen; global conf;
if ~isfield(r,'feederName'),r.feederName = 'Fallbrook';end
if isempty(r.feederName), r.feederName = fName; end
if isempty(r.feederSetup), r.feederSetup = fdSetup; end
if isempty(r.timeDay), r.timeDay = tDay; end
if isempty(r.penLevel), r.penLevel = pen; end
if ~isfield(r,'conf'), r.conf = conf; end
if ~exist(outDir,'dir'), mkdir(outDir); end
fn = sprintf('%s/%s_%s_%s_%s_Pen%03.0f',outDir,prefixStr,r.feederName,r.feederSetup,r.timeDay,r.penLevel);
end
