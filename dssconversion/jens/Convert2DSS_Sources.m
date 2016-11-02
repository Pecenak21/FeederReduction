function [cir vs] = Convert2DSS_Sources(fi,n_circuit)
% [cir vs] = Convert2DSS_Sources(fi)
%
% PURPOSE : Converts sources from object to OpenDSS
%
%
% INPUT :   cir: struct that has all entries (e.g., o.R1)
%
% OUTPUT :  one OpenDSS circuit-defining source
%           arbitrary number of voltage sources
%           Information on page 105 in OpenDSS manual V7.6(November 2012)

vs=[];

for i = 1:length(fi)
    if i==n_circuit
        cir = dsscircuit; % one circuit-defining source
        cir.Name = fi(i).NetworkID;
        cir.bus1=fi(i).NodeID;
        cir.basekv=fi(i).KVLL;
        cir.R1=fi(i).R1;
        cir.X1=fi(i).X1;
        cir.R0=fi(i).R0;
        cir.X0=fi(i).X0;
    else
        vs=dssvsource; % the rest are regular voltage sources
        vs(i).Name = ['source_' fi(i).NodeID{1}];
        vs(i).bus1=fi(i).NodeID;
        vs(i).basekv=fi(i).KVLL;
        vs(i).R1=fi(i).R1;
        vs(i).X1=fi(i).X1;
        vs(i).R0=fi(i).R0;
        vs(i).X0=fi(i).X0;
    end
end

end

