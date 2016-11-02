function [dat] = simControl(c,time,timeId,type,mode)
% simControl implements different control schemes in between consecutive simulations. 
% Control vars include PV and Energy Storge's active, reactive powers, and voltage regulator tap positions.
% Forecast info is obtained from the circuit input 'c'.
%
% Control type:
%       (1) Instantaneous 'local' control using PV's reactive powers without forecast data.
%       (2) 'local-opt' : local optimization control utilizing forecast data.
%       (3) 'distribute' : distributed optimization control utilizing forecast data.
%       (4) 'central' : centralized optimization control utilizing forecast data. 
% *Currently only support: (1).
%
% Input:
%           o   : global var that should already exists
%           c   : circuit struct
%           tStart
%           tEnd
%           type : control type. Default: 'none'
%
% Output:
%           o   : global var
%
% Example of use: 
%       dat = simControl(c,tStart,tEnd,'local')
global indent; if isempty(indent), indent = ''; end

dat = [];
if ~exist('type','var') || isempty(type)
    type = 'none'; 
end
        
switch lower(strtrim(type))
    case 'none'
        return;
    case 'local'
        if isfield(c,'pvsystem')
            % adjust reactive power based on voltage level at the interested node using the local-instantaneous control equation from Andu's PhD thesis
            vRef = 1.0; vBW = 0.01; % reference voltage and voltage bandwidth  
            for i = 1:length(c.pvsystem)
                adjustReactivePow(c.pvsystem(i),vRef,vBW);
            end
        end
    case 'local-opt' % using forecast data
        % only support 'daily' and 'yearly' mode (to obtain forecast and load data for now). If operational mode is used, the following code needs modification.
        if ismember(lower(mode),{'daily','yearly'}) 
            % loadshape id according to the pvsystem list
            [~,lspvid] = ismember([c.pvsystem.(mode)],{c.loadshape.Name}');
            % loadshape id according to the load list
            [~,lsloadid] = ismember([c.load.(mode)],{c.loadshape.Name}');
            for i = 1:length(c.pvsystem)
                % get forecast data for optimized period
                fcDat = c.loadshape(lspvid(i)).Mult(tid);
                % get load data for optimized period
                loadDat = c.loadshape(lsloadid(i)).Mult(tid);
                % run optimization
                schedule = simOpt(loadDat,fcDat,volt);
            end
        else 
            error('only support ''daily'' and ''yearly'' mode (to obtain forecast and load data for now). If operational mode is used, the following code needs modification.'); 
        end
    case 'distributed'
        algId = [2,3]; % algorithm to run
        % only works for Fallbrook feeder for now since bus0 is needed info
        % TODO: put bus0 somewhere so this would be generalized
        bus0 = '0520'; genType = 'pvsystem';
        % let pv generate reactive and storage generate active
        for i = 1:length(c.(genType))
            c.(genType)(i).pmpp = 0;                    
            c.(genType)(i).kvar = c.(genType)(i).kVA / 1.1;

            c.storage(i).kW = [];
            c.storage(i).kWrated = c.(genType)(i).kvar;
            c.storage(i).kWhrated = c.storage(i).kWrated * 2;
            c.storage(i).kvar = [];
        end
        gmat.gs = getGmatrix(c,c.storage,bus0); % get Gmatrix
        gmat.gg = getGmatrix(c,c.(genType),bus0);
        bg = getBiograph(c);
        % find neighbors of Generation/PV nodes (that also have storage systems)
        fp = fNamePrefix('Neighbor_'); fp = [fp{1} '.mat'];
        if exist(fp,'file')
            nb = load(fp);
        else
            gbus = [c.(genType).bus1]; if size(gbus,1) < size(gbus,2), gbus = cleanBus(gbus'); end
            nb.pvnb = neighbor(cleanBus(gbus), bg); % indices of neighbors based on index of c.(genType)
            nb.pvnb0 = neighbor([gbus; cleanBus(bus0)], bg); % indices of neighbors based on index of c.(genType)
            sbus = [c.(genType).bus1]; if size(sbus,1) < size(sbus,2), sbus = cleanBus(sbus'); end
            nb.ssnb = neighbor(sbus, bg);
            nb.ssnb0 = neighbor([sbus; cleanBus(bus0)], bg);
            saveFile(fp,nb); fprintf(['%sSaved neighbor nodes'' file: ' fp '\n'],indent);
        end
        % generic args
        args.plotFlag = 1; args.gs = gmat.gs; args.gg = gmat.gg; args.pvnb = nb.pvnb; args.pvnb0 = nb.pvnb0;
        args.ssnb = nb.ssnb; args.ssnb0 = nb.ssnb0; args.theta = 20; 
        args.time = time; args.timeId = timeId;
        % alg1
        fp = fNamePrefix('DistCtrlAlg1'); fp = [fp{1} '.mat'];
        if ~exist(fp,'file') && ismember(1, algId)
            niter = 10000;  args.algIndex = 1; args.saveFn = 'f520_alg1'; args.saveFreq = 10;
            d.d11 = distAlg(c,niter,'c1a1_f520',args,genType); d.args = args;
            saveFile(fp,d);
        end
        
        % alg2
        fp = fNamePrefix('DistCtrlAlg2'); fp = [fp{1} '.mat'];
        if ~exist(fp,'file') && ismember(2, algId)
            niter = 1000; args.plotFlag = 1; args.plotFreq = 10;
            args.algIndex = 2; args.gamma = 4800; args.costWeight = 1/10;
            d.d12 = distAlg(c,niter,'c1a2_f520',args,genType); d.args = args;
            saveFile(fp,d);
        end
        
        % alg3
        fp = fNamePrefix('DistCtrlAlg3'); fp = [fp{1} '.mat'];
        if ~exist(fp,'file') && ismember(3, algId)
            niter = 10000; args.plotFlag = 1; args.algIndex = 3; %args.fixNum = niter;
            d.d13 = distAlg(c,niter,'c1a3_f520',args,genType); d.args = args;
            saveFile(fp,d);
        end
    case 'centralized'
        
    otherwise
        error;
end
end

function adjustReactivePow(com,vRef,vBW,gamma)
global o;
% com: component (e.g. c.pvsystem(1))
if ~exist('vRef','var'), vRef = 1; end
if ~exist('vBW','var'), vBW = .01; end
if ~exist('gamma','var'), gamma = 0.04; end

v = dssgetval(o,com,'voltages'); % in volts
vNorm = str2double(getval(o,com,'kv'))*1000; % in volts
% voltage per unit
if com.phases == 3 
    vpu = v(1)/vNorm * sqrt(3);
else
    error('Hasn''t written code to handle component with phases different from 3 yet. FIX IT NOW');
end

if abs(vpu-vRef) > vBW/2
    % adjust reactive power based on voltage level at the interested node using the local-instantaneous control equation from Andu's PhD thesis
    pq = dssgetval(o,com,'powers'); pq = reshape(pq,2,length(pq)/2); % first row is active power, 2nd row is reactive power
    p = sum(pq(1,:))/1000;
    Qavail = sqrt( str2double(getval(o,com,'kva'))^2 - p^2 );
    q = Qavail * (1 - 2/( 1 + exp(-4/gamma*(vpu-vRef)) ));
    
    % set value in sim engine
    setval(o,com,'kvar',sprintf('%5.4f',q));
end

end

function simOpt(load,fc)
end