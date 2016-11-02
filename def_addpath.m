function def_addpath(force)
% adds the current directory and all non-.svn subdirectories to the MATLAB path
%
% this is equivalent to right-clicking in the folder browser and choosing "add
% selected folders and subfolders to path" except for we ignore .svn directories

% skip if already added the path
global PATHADDED;
if ~exist('force','var') || isempty(force) || ~force
    if PATHADDED
        return;
    end
end

%  Get the current working directory, and all its subdirectories
newpath = genpath(pwd);
%  Filter out .svn dirs
newpath = regexprep(newpath,':?[^:]*[/\\]\.svn[^:]*','');
newpath = regexprep(newpath,':?[^:]*[/\\]\.git[^:]*','');
%to the matlab path
addpath(newpath);

javaCheckInit;

PATHADDED = 1;
end