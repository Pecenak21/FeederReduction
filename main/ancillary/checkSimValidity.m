%%
simDir = normalizePath('$KLEISSLLAB24-1/database/gridIntegration/PVimpactPaperFinal');
outputDir = 'Y:/database/gridIntegration/PVimpactPaperFinal/';
x = dir([simDir '/Res*']); disp({x.name}');
x = {x.name}';
for i = 1:length(x)
    res = load([simDir '/' x{i}]); disp(x{i});
    plotSimResult(res,'single','',outputDir);
end

%% check s2b and s3b simulations
simDir = normalizePath('$KLEISSLLAB24-1/database/gridIntegration/PVimpactSimFinal');
outputDir = [simDir '/fig'];
x = dir([simDir '/Res*Fallbrook*s2b*']); disp({x.name}');
x = {x.name}';
for i = 1:length(x)
    res = load([simDir '/' x{i}]); disp(x{i});
    plotSimResult(res,'single','',outputDir);
end