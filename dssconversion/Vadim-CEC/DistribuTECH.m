function [c1] = DistribuTECH (c1, p, CaseNum)

if CaseNum==1 %changes required for case 1
    
%% change PV size
% the size in the CYME files corresponds to the momentary power flow values
% as obtained from the SCADA system. for Daily simulations we need the
% actual rated kWs.
for i=1:length(c1.pvsystem)
    c1.pvsystem(i).kVA = 500;
    c1.pvsystem(i).Pmpp = 500;
end

%% add loadshapes. taking the loadshapes from the CSI project. 
% these loadshapes are used for the DistribuTECH Brazil paper. for the actual
% CEC results different loadshapes will be used.
loadshape_resident = dssparse('LoadShape-Residential(from CSI).dss');
loadshape_resident = loadshape_resident.loadshape;
loadshape_small = dssparse('LoadShape-Small(from CSI).dss');
loadshape_small = loadshape_small.loadshape;
loadshape_medium = dssparse('LoadShape-Medium(from CSI).dss');
loadshape_medium = loadshape_medium.loadshape;
loadshape_large = dssparse('LoadShape-Large(from CSI).dss');
loadshape_large = loadshape_large.loadshape;
loadshape_PV = dssparse('PVShape(from CSI).dss');
loadshape_PV = loadshape_PV.loadshape;
loadshape=[loadshape_resident, loadshape_small, loadshape_medium, loadshape_large, loadshape_PV];
[c1.loadshape]=deal(loadshape);

%% increase the loadshapes
for ii=1:length(c1.loadshape)
    AA=c1.loadshape(ii).Mult;
    Max=max(AA);
    c1.loadshape(ii).Mult=AA./Max;
end

%% additional changes
for vz=1:length(c1.load)
    c1.load(vz).Daily = ['loadshape_' c1.load(vz).Name];  
end
for zz=1:length(c1.pvsystem)
    c1.pvsystem(zz).daily = ['pv' num2str(zz)];
%     c1.pvsystem(zz).kv = 12;% change pv voltage from 0.208 to 12kv (if decided to remove trafos)
end
% for zz=1:length(c1.transformer)-1
%     c1.transformer(zz).kVs = [12 12];
% end
while length(c1.regcontrol) > 1
    c1.regcontrol(1) = [];
end

%% add voltage regulators (to provide a voltage control effect)
% [vz,vz] = ismember('105841142_03265_GS1509_4_03265', c1.line(:).Name);
% n = c1.transformer(end);
% n.Name = 'reg1';
% n.Conns = {'wye' 'wye'};
% n.Buses = {c1.line(vz).bus1 c1.line(vz).bus2};
% n.kVas = [20000 20000];
% n.kVs = [12 12];
% c1.transformer(end+1) = n;
% c1.line(vz) = '';
% n = c1.regcontrol(end);
% n.Name = c1.transformer(end).Name;
% n.transformer = c1.transformer(end).Name;
% n.vreg = 100;
% n.band = 3;
% n.winding = 2;
% c1.regcontrol(end+1) = n;
% c1.transformer(end).XHL = '';
% c1.transformer(end).sub = '';
% c1.transformer(end).Rs = [];
% 
% %% changes to substation Trafo
% c1.regcontrol(1).vreg = 120;

%% add monitors
n = dssmonitor('Name','Sub');% high voltage side
n.Element = 'Transformer.SAN_BERDO';
n.Mode = 0;
% Element = 1;
c1.monitor = n;

for yy=1:length(c1.line)
    n = dssmonitor('Name',[c1.line(yy).Name '.' c1.line(yy).bus1]);
    n.Element = ['Line.' c1.line(yy).Name];
    n.Mode = 0;
    n.Terminal = 1;
    c1.monitor(end+1) = n;
end

% SAVE MONITORS - OPTIONAL if you want to save all monitors as csv files
% after performing daily simulaitons. The DailySimulaiton.m runs the
% SCE_Centaur.dss file to export all monitors.
monitor_dir = [p(1:strfind(p,'/SCE_Centaur.dss')) 'Monitors'];
if ~exist([p(1:strfind(p,'/SCE_Centaur.dss')) 'Monitors'], 'dir');
    monitor_dir = [p(1:strfind(p,'/SCE_Centaur.dss')) 'Monitors'];
    mkdir(monitor_dir);
end

dfn = 'Centaur_mon.dss';
for BB=2:length(c1.monitor)
    fid = fopen([monitor_dir '/' dfn], 'at');
    if(fid==-1), error('dsswrite:openfailed','Failed to open output file %s for writing!\nRemember to close open files before overwriting.',[savepath '/' dfn]); end
    str = ['export monitor ' c1.monitor(BB).Name];
    try
        fprintf(fid,'%s\n', str);
        fclose(fid);
    catch err
        warning('dsswrite:openfiles','Remember to close files before overwriting them!');
        rethrow(err);
    end
end
else %changes required for case 2
%     [vz,vz] = ismember('loadshape_lP5506157_T_03265', c1.loadshape(:).Name);
% %     c1.loadshape(vz).Mult(1200:1320) = zeros(1,121);
%     c1.loadshape(vz).Mult(1200:1320) = zeros(1,121);
    [vz,vz] = ismember('loadshape_lP5517477_T_03265', c1.loadshape(:).Name);
    c1.loadshape(vz).Mult(1200:1320) = zeros(1,121);
    [vz,vz] = ismember('loadshape_lP5517478_T_03265', c1.loadshape(:).Name);
    c1.loadshape(vz).Mult(1200:1320) = zeros(1,121);    
end
end