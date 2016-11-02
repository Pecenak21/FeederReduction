function [ data ] = dssSimulationTMP( circuit, mode, t, doPlotVoltProf)
%DSSSIMULATION Run simulation and store data in circuitout.data
%   [ circuitout ] = dssSimulation( circuit, mode, steptime)
%   INPUTS:
%       "circuit": circuit on which run simmulation             [struct]
%       "mode": simulation mode                                 [char]
%       "steptime": steptime in SECONDS for daily and yearly modes   [double]
%   OUTPUT:
%       "circuitout": circuit with data field
%   [ circuitout ] = dssSimulation( ct, 'daily', 30)

%% use loadshape for capacitors for simulations with PV
global indent; if isempty(indent), indent = ''; end
global conf;
if isfield(circuit, 'dataNoPV')
    capShape = 1;
else
    capShape = 0;
end

if ~exist('conf','var') || isempty(conf)
    conf.ctrlType = 'none';
    conf.ctrlHorizon = 15; % in minutes
    conf.ctrlTStep = 30; % in seconds
else
    x = conf;
end

%% Run Simulation
% Set simulation
if ~exist('mode','var') || isempty(mode)
    mode = 'snapshot';
end

% add energy meter
circuit = addEnergyMeter(circuit);
% p = dsswrite(circuit,[],0,[]);
p=circuit;
% run simulation
global o;
o = actxserver('OpendssEngine.dss');
o.Start(0);
dssText = o.Text;
global oStartTime;
oStartTime = t(1); nt = length(t);

%% OpendDSS simulation with No PV included
dssText.Command = 'Clear';
cDir = pwd;
dssText.Command = ['Compile "' p '"'];
cd(cDir);
dssText.Command = 'Set controlmode = static';
dssText.Command = ['Set mode = ' mode ];
dssText.Command = ['Set stepsize = ' num2str(conf.timeStep) 's'];
dssText.Command = 'Set number = 1'; % number of time steps or solutions to run or the number of Monte Carlo cases to run.
dssText.Command = 'SetLoadandGenkV';

fprintf('%sRunning simulation...\n',indent)
dssCircuit = o.ActiveCircuit;
Dist = dssCircuit.AllNodeDistances; % Distance in km
nodeName = dssCircuit.AllNodeNames;
% dssEnergyMeters = dssCircuit.Meters;
% Dist = dssCircuit.AllNodeDistances*(3280.8399); % conversion from km to kft (not sure why OpenDSS presents it in km, when the line lengths are clearly identified to be in kft)

dssSolution = dssCircuit.Solution;
dssSolution.MaxControlIterations=100;
dssSolution.MaxIterations=100;
dssSolution.InitSnap; % Initialize Snapshot solution
dssSolution.dblHour = 0.0;

% initialize output
Volt_MaxMin=zeros(nt,2); LossTotal=zeros(nt,2); LossLine=zeros(nt,2);
Volt=zeros(nt,length(dssCircuit.AllBusVmagPu)); TotalPower=zeros(nt,2);
% Volt_NoPV_PhaseA = []; Volt_NoPV_PhaseB = []; Volt_NoPV_PhaseC = [];
converged = zeros(nt,1); tapPos = struct;
if isfield(circuit,'pvsystem')
    PVpwr.kw = zeros(nt,length(circuit.pvsystem));
    PVpwr.kvar = zeros(nt,length(circuit.pvsystem));
end

simHours = (t(end) - t(1))*24; % simulation time in hours
tSim = t-t(1); % simulation time;
while dssSolution.dblHour <= simHours
    [~,i] = ismember(dssSolution.dblHour,tSim*24);
    if isempty(i)
        error('Simulation time is not in the desired array of time to run! Check. Possible issue: Rounding error!');
    end
    if i > 1 % apply control from the second time step
        if exist('conf','var') && isfield(conf,'ctrlType') && ~isempty(conf.ctrlType)
            % get time index of interest using forecast info for control
            tend = getTimeId(t(i)+conf.fcTimeAhead/60/24,t);
            if isempty(tend), error('Time not in the specified array!'); end;
            dat = simControl(circuit,t(i:tend),i:tend,conf.ctrlType,mode);
        end
    end
    if mod(dssSolution.dblHour*3600,3600)==0
        fprintf('%sId: %04d Time: %02d h %04.0f s ... \n',indent,i,floor(dssSolution.dblHour),mod(dssSolution.dblHour*3600,3600));
    end
    if capShape
        idx = round(dssSolution.dblHour*3600/circuit.loadshape(2).sInterval)+1;
        for j = 1:length(circuit.capacitor),
            kvar = circuit.capacitor(j).kvar*circuit.loadshape(2).Mult(idx);
            dssText.Command = ['Capacitor.' circuit.capacitor(j).name '.kvar = ' num2str(kvar)];
        end
    end
    dssSolution.Solve;
    if dssSolution.Converged
        converged(i) = 1;
        if i+1 <= length(tSim), o = setSolutionTime(o,t(i+1)); end
        % record output when converged
        Volt_MaxMin(i,:) = [max(dssCircuit.AllBusVmagPu(dssCircuit.AllBusVmagPu < 4)) min(dssCircuit.AllBusVmagPu(dssCircuit.AllBusVmagPu > 0.3))];
        LossTotal(i,:) = dssCircuit.Losses/(1e+03); % MW/MVar
        Volt(i,:) = dssCircuit.AllBusVmagPu;
        if exist('doPlotVoltProf','var') && doPlotVoltProf
            plotVoltProf(Volt(i,:),Dist,t(i),i);
        end
        LossLine(i,:) = dssCircuit.LineLosses;
        TotalPower(i,:) = dssCircuit.TotalPower/(1e+03)*(-1);
        %get power production of pvsystem
            if isfield(circuit, 'pvsystem')
                pw = dssgetval(o,circuit.pvsystem,'powers');
                pvV = dssgetval(o,circuit.pvsystem,'voltages');
				PVpwr.kw(i,:) = (pw(:,1)+pw(:,3)+pw(:,5)+pw(:,7))';
                PVpwr.kvar(i,:) = (pw(:,2)+pw(:,4)+pw(:,6)+pw(:,8))';
				PVpwr.V(i,:) = (pvV(:,1)+pvV(:,3)+pvV(:,5)+pvV(:,7))';
			end
		% get voltage regulator tap positions from the transformer
        val = getval(o,circuit.transformer,'taps');
        tPos = readTranxTap(val);
        for j = 1:length(circuit.transformer)
            tapPos.transformer(j).V(i,:) = tPos(j,:);
        end
        % get capacitor
        if isfield(circuit,'capacitor')
            cPos = readCapTap(getval(o,circuit.capacitor,'states'));
            for j = 1:length(circuit.capacitor)
                tapPos.capacitor(j).pos(i,:) = int32(cPos{j});
            end
        end
    else
        fprintf('%sId: %04d Time: %d h %4.0f s ... ',indent,i,floor(dssSolution.dblHour),mod(dssSolution.dblHour*3600,3600));
        setSolutionTime(o,t(i));
        if dssSolution.MostIterationsDone < dssSolution.MaxIterations &&  dssSolution.ControlIterations == dssSolution.MaxControlIterations && ~dssSolution.ControlActionsDone
            warning('%sWARNING: System converged but control devices did not converge (usually fine still). Will continue.\n',indent);
            converged(i) = 2;
        else
            converged(i) = 0;
            error('%sALARM: Solution doesn''t converge! Trying to continue anyway! No data is recorded! Check your circuit again! Will try to continue.\n', indent);
        end
    end
    
end

if isfield(circuit,'transformer')
    % post process tap position
    maxTap = str2double(getval(o,circuit.transformer,'maxtap'));
    minTap = str2double(getval(o,circuit.transformer,'mintap'));
    numTaps = str2double(getval(o,circuit.transformer,'numtaps'));
    % calculate transformer one tap's voltage change
    dTap = (maxTap-minTap)./numTaps;
    for i = 1:length(circuit.transformer)
        tapPos.transformer(i).dTap = dTap(i);
        % calculate position from voltage levels
        tapPos.transformer(i).pos = (tapPos.transformer(i).V-1)/dTap(i);
        % calculate number of tap changes during period of simulation
        tapPos.transformer(i).numTapChange = sum(abs(diff(tapPos.transformer(i).pos)));
    end
    data.totTranxTapOpe = nansum([tapPos.transformer.numTapChange]);
else
    data.totTranxTapOpe = 0;
end

% process cap controller data
if isfield(circuit,'capacitor')
    cSize = cellfun(@sum,readCapTap(getval(o,circuit.capacitor,'kvar')));
    if length(circuit.capacitor) == 1
        nStep = str2double(getval(o,circuit.capacitor,'numSteps'));
    else
        nStep = cellfun(@str2double,getval(o,circuit.capacitor,'numSteps'));
    end
    for i = 1:length(circuit.capacitor)
        tapPos.capacitor(i).sumSubCap = sum(tapPos.capacitor(i).pos,2);
        tapPos.capacitor(i).numTapChange = sum(sum(abs(diff(tapPos.capacitor(i).pos))));
        tapPos.capacitor(i).numStep = nStep(i);
        tapPos.capacitor(i).capSizeKvar = cSize(i);
    end
    data.totCapTapOpe = sum([tapPos.capacitor.numTapChange]);
else
    data.totCapTapOpe = 0;
end

data.nodeName = nodeName';
data.Voltage = Volt;
data.LineLoss = LossLine;
data.TotalLoss = LossTotal;
data.VoltMaxMin = Volt_MaxMin;
data.TotalPower = TotalPower;
data.Dist=Dist;
data.tapPos = tapPos;
data.converged = converged;
data.numConverged = sum(converged==1);
data.numNotConverged = sum(converged==0);
data.numMaxCtrlExceed = sum(converged==2);
data.Pvdata=PVpwr;
data.time = t;
% if isfield(circuit,'pvsystem')
%     data.PVpwr = PVpwr;
% end

delete(o);
end

function fp = plotVoltProf(v,dist,t,i)
figure(100);
plot(dist,v,'.');
xlabel('Distance from substation [km]');
ylabel('Voltage [p.u.]');
title(datestr(t,'yyyy-mm-dd HH:MM:SS'));
grid on; box on;
fp = ['tmp/voltProf' '_' sprintf('%04d',i) '.png'];
saveas(gcf,fp);
fprintf('Saved file: %s\n',fp);
end