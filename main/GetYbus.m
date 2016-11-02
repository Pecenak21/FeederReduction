function [Ybus, YbusOrderVect, YbusPhaseVect, buslist] = GetYbus(circuit)

%Created by Zachary K. Pecenak on 6/18/2016

%Example input
% load('c:\users\zactus\gridIntegration\results\ValleyCenter_wpv_existing.mat')
% Bus_upstream={'03551325'};
% or Bus_upstream=c.buslist.id(ceil(length(c.buslist.id)*rand(1,1))) this is a random bus
% [c] = FeederReduction(Bus_upstream,c);
%Check both inputs are met

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% get Ybus from OpenDSS for feeder.
fprintf('\nGetting YBUS and buslist from OpenDSS: ')
tic

%load the circuit and generate the YBUS
p = dsswrite(circuit_woPV,[],0,[]); o = actxserver('OpendssEngine.dss');
dssText = o.Text; dssText.Command = 'Clear'; cDir = pwd;
dssText.Command = ['Compile "' p '"']; dssCircuit = o.ActiveCircuit;
Ybus=dssCircuit.SystemY;

%Convert the Ybus to a matrix
ineven=2:2:length(Ybus); inodd=1:2:length(Ybus);  
Ybus=Ybus(inodd)+1i*Ybus(ineven); Ybus=reshape(Ybus,sqrt(length(Ybus)),sqrt(length(Ybus)));
Ybus=sparse(Ybus);
%get buslist in order of Ybus and rearrange
busnames=regexp(dssCircuit.YNodeOrder,'\.','split');
YbusOrderVect=[busnames{:}]'; YbusOrderVect(find(cellfun('length',YbusOrderVect)==1))=[];
YbusPhaseVect=[busnames{:}]'; YbusPhaseVect(find(cellfun('length',YbusPhaseVect)>1))=[]; YbusPhaseVect=str2double(YbusPhaseVect);
% %Here we generate the list of node numbers in the order of the Ybus
buslist=dssCircuit.AllBusNames;
clear inodd ineven 

t_=toc;
fprintf('time elapsed %f',t_)
end