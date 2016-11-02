function [Start,Obj,Text] = dssStartup
% Function for starting up the dss

%instantiate the dss Object
Obj = actxserver('OpendssEngine.dss');
%
%Start the dss.   Only needs to be executed the first time w/in a
%Matlab session
Start = Obj.Start(0);

% Define the text interface
Text = Obj.Text; 
end