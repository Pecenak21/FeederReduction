function [ d ] = distAlg( c, niter, name, args, genType )
% Control simulation algorithm option (args.algIndex):
%
% (1) No hard forcing on contraints of Gen & SS (No cvx toolbox needed).
% Let the algorithm converges naturally. Takes long time but will get to
% final solution.
%
% (2) (Default) Use CVX optimization toolbox to force constraints of generations of kvar (PV systems) and kw (storage
% systems) before scheduling. No Model Predictive Control used. Using Optical Power Flow to find
% solution for power flow equations and schedule components from then.
%
% (3) Use CVX optimization toolbox to force constraints of generations of kvar (PV systems) and kw (storage
% systems) before scheduling. Model Predictive Control is used to replace
% power flow equations -> True distributed algorithm.
%
% (4)
%
% args.gs fields:  algIndex (algorithm index [1,2,3], default:2)
%               args.gg (gen matrix), args.gs (storage matrix), bg (biograph),
%               args.pvnb (pv neighbor), args.pvnb0 (with bus0)
%               args.ssnb (storage neighbor), args.ssnb0 (with bus0)
%               fixIterNum (fixed number of iterations)
%               costWeight (weight of power cost relative to lineloss), Default: 1/U0

% create a defacto circuit just for this algorithm since we're not simulation in time series but try to find a solution using dssget
% loadshapes are not needed
% ct = rmfield(c,'loadshape'); 

global conf;

if ~exist('niter','var') || isempty(niter), niter = 10000; end
try args.algIndex; catch, args.algIndex = 2; end
if ~exist('name','var') || isempty(name), name = ['alg' num2str(args.algIndex)]; end

% Theta
try theta = args.theta; catch
    try theta = 180/pi*atan2( c.linecode.Xmatrix, c.linecode.Rmatrix);
    catch, theta = 20; end
end
try Wmax = args.Wmax; catch, Wmax = 1.05^2; end
try Wmin = args.Wmin; catch, Wmin = 0.95^2; end
try args.fixIterNum; catch, args.fixIterNum = 0; end

% stats
% gamma value -> large means fast convergence but watch out for
% non-convergence when gamma is too large
try y = args.gamma; catch, y = 10^-1; end
% storage initial condition
try x0 = args.x0; catch, x0 = 1; end
try powerCost = args.powerCost; catch 
    powerCost = ones(length(args.time),1); 
end
%tolerance
try tor = args.tor; catch, tor = 0.01/100; end % in percents
% cost weight
U0 = c.circuit.basekv * 1e3;
try costWeight = args.costWeight; catch, costWeight = 1/U0; end

% %% initialize bat charges from simulation results
try args.time; catch, args.time = 1:5; end
nt = length(args.time);
nS = length(c.storage);
nG = length(c.(genType));

d.costtotal = nan(niter,1);
d.subPower = complex(zeros(niter,nt),0);
d.subLoss = complex(zeros(niter,nt),0);
if isfield(d,'nuvec') && niter > size(d.nuvec,1)
    xS = zeros(niter-size(d.nuvec,1),nS,nt);
    xG = zeros(niter-size(d.nuvec,1),nG,nt);
    xS1 = zeros(niter-size(d.nuvec,1),nS+1,nt);
    xG1 = zeros(niter-size(d.nuvec,1),nG+1,nt);
    
    d.gluvec = [d.gluvec; xG];
    d.gllvec = [d.gllvec; xG];
    
    d.qvec = [d.qvec; xG];
    d.vvec = [d.vvec; xS];
    
    d.usvec = [d.usvec; xS1];
    d.asvec = [d.asvec; xS1];
    d.ugvec = [d.ugvec; xG1];
    d.agvec = [d.agvec; xG1];
    d.wgvec = [d.wgvec; xG1];
    
    switch args.algIndex
        case 1
            d.nuvec = [d.nuvec; xS];
            d.nlvec = [d.nlvec; xS];
            d.muvec = [d.muvec; xS];
            d.mlvec = [d.mlvec; xS];
        otherwise
    end
    kStart = size(d.nuvec,1);
else
    kStart = 2;
    d.gluvec = zeros(niter,nG,nt); % lamdba
    d.gllvec = zeros(niter,nG,nt); % lamdba
    
    d.qvec = zeros(niter,nG,nt);
    d.vvec = zeros(niter,nS,nt);
    
    % storage systems' voltages with node 0 (bus 799)
    d.usvec = zeros(niter,nS+1,nt);
    d.asvec = zeros(niter,nS+1,nt);
    % %% initialize gen voltages without node 0
    d.ugvec = zeros(niter,nG+1,nt);
    d.agvec = zeros(niter,nG+1,nt);
    d.wgvec = zeros(niter,nG+1,nt);
    % 2D data
    
    glu = zeros(nG,nt); % lamdba
    gll = zeros(nG,nt); % lamdba
    
    q = zeros(nG,nt);
    v = zeros(nS,nt);
    
    switch args.algIndex
        case 1
            d.nuvec = zeros(niter,nS,nt); % upper bound of nl
            d.nlvec = zeros(niter,nS,nt); % nl lower bound
            d.muvec = zeros(niter,nS,nt);
            d.mlvec = zeros(niter,nS,nt);
            
            nu = zeros(nS,nt); % upper bound of nl
            nl = zeros(nS,nt); % nl lower bound
            mu = zeros(nS,nt);
            ml = zeros(nS,nt);
            
            % 2D data previous (old)
            nu_o = zeros(nS,nt); % upper bound of nl
            nl_o = zeros(nS,nt); % nl lower bound
            mu_o = zeros(nS,nt);
            ml_o = zeros(nS,nt);
        otherwise
    end
    glu_o = zeros(nG,nt); % lamdba
    gll_o = zeros(nG,nt); % lamdba
    
    q_o = zeros(nG,nt);
    v_o = zeros(nS,nt);
    
    % loadshape dat
    lsdat = [c.loadshape.Mult]';
    % loadshape ids for generators/ pvsystems
    [~,lsidgen] = ismember([c.(genType).(conf.mode)],{c.loadshape.Name}');
    % generators' active power
    genP = lsdat(lsidgen,:);
    
    % loadshape ids for loads
    [~,lsidload] = ismember([c.load.(conf.mode)],{c.loadshape.Name}');
    loadPQ = lsdat(lsidload,:);
    
    % first simulation to initialize some variables
    [vol,~,da] = controlSim(c,genType,v,q,args.time,[],[] );%%
    us_o = vol.us; ug_o = vol.ug; as_o = vol.as; ag_o = vol.ag; wg_o = vol.wg;
%     us(1,:,:) = vol.us;
%     ug(1,:,:) = vol.ug;
%     as(1,:,:) = vol.as;
%     ag(1,:,:) = vol.ag;
%     wg(1,:,:) = vol.wg;
    
    d.costtotal(1) = calCost(da,costWeight);
    d.subPower(1,:) = da.subPower';
    d.subLoss(1,:) = da.subLoss';
end

switch(args.algIndex)
    case {2,3}
        cvx_startup; cvx_setup;
    otherwise
end
%% SIMULATION
tx = tic;
switch(args.algIndex)
    case 3
        lmaxvec = zeros(1,niter);
    otherwise
end
for k = kStart:niter
    tic; disp(k);
    %% calculate coefficients
    for t = 1:nt
        %% generation network in which each generator is accompanied by a storage system
        for l = 1:nG
            glu(l,t) = max(0, glu_o(l,t) + y * (wg_o(l,t) - Wmax) );
            gll(l,t) = max(0, gll_o(l,t) + y * (Wmin - wg_o(l,t)) );
            
            switch args.algIndex
                case 1
                    nu(l,t) = max(0, nu_o(l,t) + y * (v_o(l,t) - (c.storage(l).kWrated*1000) )/(2*(c.storage(l).kWrated*1000)) );
                    nl(l,t) = max(0, nl_o(l,t) + y * (- (c.storage(l).kWrated*1000) - v_o(l,t) )/(2*(c.storage(l).kWrated*1000)) );
                    mu(l,t) = max(0, mu_o(l,t) + y * ( x0 + 1/(c.storage(l).kWhrated*1000) * sum(reshape(v_o(l,1:t),1,t)) - 1 ) );
                    ml(l,t) = max(0, ml_o(l,t) + y * ( -x0 - 1/(c.storage(l).kWhrated*1000) * sum(reshape(v_o(l,1:t),1,t)) ) );
                otherwise
            end
        end
    end
    
    %% calculate active charging rate for storage and reactive generation for generators
    for t = 1:nt
        for l = 1:nG
            gid = find(args.pvnb.neighborMatrix(l,:));
            gid0 = find(args.pvnb0.neighborMatrix(l,:));
            % id of neigbors
            nid = find(args.ssnb.neighborMatrix(l,:));
            nid0 = find(args.ssnb0.neighborMatrix(l,:));
            
            A = 0;
            for z = 1:length(nid0)
                h = nid0(z);
                A = A + args.gs(l,h)* us_o(l,t) * us_o(h,t) * cosd( as_o(h,t) - as_o(l,t) - theta );
            end
            
            B = 0;
            for z = 1:length(nid)
                h = nid(z);
                switch args.algIndex
                    case 1
                        B = B + args.gs(l,h)* ( -(nu(h,t) -nl(h,t))/(2*(c.storage(l).kWrated*1000)) - 1/(c.storage(l).kWhrated*1000) * ...
                            sum( reshape( mu(h,:),1,nt ) - reshape( ml(h,:),1,nt ) ) ...
                            - costWeight*powerCost(t) );
                    otherwise
                        B = B + args.gs(l,h)* ( - costWeight*powerCost(t) );
                end
                
            end
            v(l,t) = v_o(l,t) + cosd(theta)*(glu(l,t)-gll(l,t)) + A + U0^2/2 * B;
            
            
            C = sum( args.gg(l,gid0) .* ug_o(l,t) .* ug_o(gid0,t)' .* sind( ag_o(gid0,t) - ag_o(l,t) -theta)');
            D = sum(args.gg(l,gid));
            q(l,t) = q_o(l,t) - sind(theta) * (glu(l,t) - gll(l,t) ) + C + U0^2/2 * costWeight* powerCost(t) * sind(theta) * D;
        end
        
        for l = nG+1:nS
            % id of neigbors
            nid = find(args.ssnb.neighborMatrix(l,:));
            nid0 = find(args.ssnb0.neighborMatrix(l,:));
            
            A = 0;
            for z = 1:length(nid0)
                h = nid0(z);
                A = A + args.gs(l,h)* us_o(l,t) * us_o(h,t) * cosd( as_o(h,t) - as_o(l,t) - theta );
            end
            %%
            B = 0;
            for z = 1:length(nid)
                h = nid(z);
                switch args.algIndex
                    case 1
                        B = B + args.gs(l,h)* ( -(nu(h,t) -nl(h,t))/(2*(c.storage(l).kWrated*1000)) - 1/(c.storage(l).kWhrated*1000) * ...
                            sum( reshape( mu(h,:),1,nt ) - reshape( ml(h,:),1,nt ) ) ...
                            - costWeight*powerCost(t) );
                    otherwise
                        B = B + args.gs(l,h)* ( - costWeight*powerCost(t) );
                end
            end
            v(l,t) = v_o(l,t) + A + U0^2/2 * B;
        end
    end
    if any(isnan(v(:))) || any(isinf(v(:)))
        error('ERROR! The algorithm blew up. Check gamma and power weight values again! Suargs.ggestion: reduce them!');
    end
    
    switch args.algIndex
        case 1
        case {2,3}
            % apply projection to find viable charging state for each storage system
            for i = 1:nS
                v(i,:) = optBat(v(i,:)',x0,[-c.storage(i).kWrated, c.storage(i).kWrated]*10^3,args.time,c.storage(i).kWhrated*10^3);
            end
            % apply projection find viable reactive generation state
            for i = 1:nG
                %q(i,:) = optGen(q(i,:)', repmat(c.(genType)(i).Kw*1000,nt,1), repmat(c.(genType)(i).kVA*1000,nt,1),nt);
                q(i,:) = optGen( q(i,:)', genP(i,args.timeId)', repmat(c.(genType)(i).kVA/1.1*1000,nt,1) , nt );
            end
        otherwise
    end
    %% rerun simulation to find u(:,:)
    switch args.algIndex
        case {1,2}
            [vol, ~, da] = controlSim(c,genType,v,q,args.time,[],[]);
        case 3
            %% MPC calculation
            lmax = 10000;
            delRe = zeros(nS,nt,lmax); delIm = delRe;
            betaRe = zeros(nS,nt); betaIm = betaRe;
            u = zeros(nS,nt);
            
            if ~exist('u0','var')
                u0 = us_o .* ( cosd(as_o) +1i*sind(as_o) );
            end
            s0 = nan(1,nS); s = nan(nS,nt);
            for l = 1:nS
                if l <= nG % l is a G node
                    s0(l) = 0;
                else
                    % storage node (right now there is no loads attached so p(l,t) = 0)
                    s0(l) = 0;
                end
                
                for t = 1:nt
                    delRe(l,t,1) = 0;
                    delIm(l,t,1) = 0;
                    
                    if l <= nG % l is a G node
                        s(l,t) = 0 - v(l,t) + 1i*q(l,t);
                    else
                        % storage node (right now there is no loads attached so p(l,t) = 0)
                        s(l,t) = 0 - v(l,t) + 0;
                    end
                    
                    % voltage angle at sourcebus (701) is 0
                    phi = 0;
                    betaRe(l,t) = cosd(theta + phi)/U0 * real(conj(s(l,t)) - conj(s0(l))) - sind(theta + phi)/U0 * imag(conj(s(l,t)) - conj(s0(l)));
                    betaIm(l,t) = sind(theta + phi)/U0 * real(conj(s(l,t)) - conj(s0(l))) + cosd(theta + phi)/U0 * imag(conj(s(l,t)) - conj(s0(l)));
                end
                
                % h < 2*|G U S| => h < 2*9=18
                h = 1;
                % id of neigbors without node 0
                nid = setdiff(find(args.ssnb.neighborMatrix(l,:)),l); % excluding itself and node 0
                % 		nid0 = find(args.ssnb0.neighborMatrix(l,:));
                for m = 2:lmax
                    for t=1:nt
                        delRe(l,t,m) = (1-h)*delRe(l,t,m-1) - h/args.gs(l,l) * ( sum( args.gs(l,nid) .* delRe(nid,t,m-1)' ) - betaRe(l,t) );
                        delIm(l,t,m) = (1-h)*delIm(l,t,m-1) - h/args.gs(l,l) * ( sum( args.gs(l,nid) .* delIm(nid,t,m-1)' ) - betaIm(l,t) );
                    end
                    if mod(m,1000)==0
                        figure(102), plot(reshape(delRe(1,1,1:m),1,m));
                    end
                    if max(max(abs( delRe(:,:,m) - delRe(:,:,m-1) ))) < tor && max(max(abs( delIm(:,:,m) - delIm(:,:,m-1) ))) < tor
                        lmax = m;
                        figure(102), plot(reshape(delRe(1,1,1:m),1,m));
                        break;
                    end
                end
                lmaxvec(k) = lmax;
                
                for t = 1:nt
                    u(l,t) = delRe(l,t,lmax) + 1i*delIm(l,t,lmax) + u0(l,t);
                end
            end
            
            %                 vol_OPF = vol;
            vol.us(1:end-1,:) = abs(u); vol.as(1:end-1,:) = angle(u)*180/pi;
            vol.ug(1:end-1,:) = vol.us(1:nG,:); vol.ag(1:end-1,:) = vol.as(1:nG,:);
            vol.wg(1:end-1,:) = vol.ug(1:end-1,:).^2./U0^2;
            
            %% TODO, how to calculate power loss and consumption for MPC case
        otherwise
    end
    
    d.costtotal(k) = calCost(da,costWeight);
    d.subPower(k,:) = da.subPower';
    d.subLoss(k,:) = da.subLoss';
    
    % update voltage profiles
    d.usvec(k,:,:) = vol.us;
    d.asvec(k,:,:) = vol.as;
    d.ugvec(k,:,:) = vol.ug;
    d.agvec(k,:,:) = vol.ag;
    d.wgvec(k,:,:) = vol.wg;
    d.vvec(k,:,:) = v;
    d.qvec(k,:,:) = q;
    
    switch args.algIndex
        case 1
            d.nuvec(k,:,:) = nu;
            d.nlvec(k,:,:) = nl;
            d.muvec(k,:,:) = mu;
            d.mlvec(k,:,:) = ml;
        otherwise
    end
    
    d.gluvec(k,:,:) = glu;
    d.gllvec(k,:,:) = gll;
    
    plotting(args,k,d);
%     saving(args,d,niter,k)
    toc;
    
    % check converging conditions
    if ~args.fixIterNum
        if max(max(abs(v-v_o)))/U0 < tor && max(max(abs(q-q_o)))/10^3 < tor
            converge = 1;
            break;
        end
    end
    us_o = vol.us; ug_o = vol.ug; as_o = vol.as; ag_o = vol.ag; wg_o = vol.wg;
    v_o = v; q_o = q; glu_o = glu; gll_o = gll;
    switch args.algIndex
        case 1
            nu_o = nu; nl_o = nl; mu_o = mu; ml_o = ml;
        otherwise
    end
end
toc(tx);

switch args.algIndex
    case 1
        d.nuvec = nuvec;
        d.nlvec = nlvec;
        d.muvec = muvec;
        d.mlvec = mlvec;
    otherwise
end

d.k = k;
d.costWeight = costWeight; % weight of cost of power
d.costPower = powerCost; % cost of power
d.niter = niter;
d.nt = nt;
d.tor = tor;
try d.converge = converge;
catch 
    if args.fixIterNum
        d.converge = nan;
    else
        d.converge = 0;
    end
end
d = cleanOutput(d);
end

function plotting(args,k,d)
if ~args.plotFlag || mod(k,args.plotFreq) > 0
    return;
end
switch args.algIndex
    case 1
        figure(100),title('All Generators');
        figure(100), subplot(3,3,1); plot(d.nuvec(1:k,:,1)); title('nu, S=1, t=1');
        figure(100), subplot(3,3,2); plot(d.nlvec(1:k,:,1)); title('nl, S=1, t=1');
        figure(100), subplot(3,3,3); plot(d.muvec(1:k,:,1));title('mu, S=1, t=1');
        figure(100), subplot(3,3,4); plot(d.mlvec(1:k,:,1));title('ml, S=1, t=1');
        figure(100), subplot(3,3,5); plot(d.gluvec(1:k,:,1));title('\lambda_U, S=1, t=1');
        figure(100), subplot(3,3,6); plot(d.gllvec(1:k,:,1));title('\lambda_L, S=1, t=1');
        figure(100), subplot(3,3,7); plot(d.qvec(1:k,:,1)/1000);title('Reactive support: q, S=1, t=1'); ylabel('kVar');
        figure(100), subplot(3,3,8); plot(d.vvec(1:k,:,1)/1000);title('Charging rate: v, S=1, t=1'); ylabel('kW');
        figure(100), subplot(3,3,9); plot(d.ugvec(1:k,:,1)/U0);title('U, G=1, t=1'); ylabel('Voltage, p.u.');
        
        figure(101), title('All Times');
        figure(101), subplot(3,3,1); plot(reshape(d.nuvec(1:k,1,:),k,nt)); title('nu, S=1, t=1');
        figure(101), subplot(3,3,2); plot(reshape(d.nlvec(1:k,1,:),k,nt)); title('nl, S=1, t=1');
        figure(101), subplot(3,3,3); plot(reshape(d.muvec(1:k,1,:),k,nt));title('mu, S=1, t=1');
        figure(101), subplot(3,3,4); plot(reshape(d.mlvec(1:k,1,:),k,nt));title('ml, S=1, t=1');
        figure(101), subplot(3,3,5); plot(reshape(d.gluvec(1:k,1,:),k,nt));title('\lambda_U, S=1, t=1');
        figure(101), subplot(3,3,6); plot(reshape(d.gllvec(1:k,1,:),k,nt));title('\lambda_L, S=1, t=1');
        figure(101), subplot(3,3,7); plot(reshape(d.qvec(1:k,1,:),k,nt)/1000);title('Reactive support: q, S=1, t=1'); ylabel('kVar');
        figure(101), subplot(3,3,8); plot(reshape(d.vvec(1:k,1,:),k,nt)/1000);title('Charging rate: v, S=1, t=1'); ylabel('kW');
        figure(101), subplot(3,3,9); plot(reshape(d.ugvec(1:k,1,:),k,nt)/U0);title('U, G=1, t=1'); ylabel('Voltage, p.u.');
    case {2,3}
        figure(100),title('All Generators');
        subplot(2,3,1); plot(d.gluvec(1:k,:,1));title('\lambda_U, S=1, t=1');
        subplot(2,3,2); plot(d.gllvec(1:k,:,1));title('\lambda_L, S=1, t=1');
        subplot(2,3,4); plot(d.qvec(1:k,:,1)/1000);title('Reactive support: q, S=1, t=1'); ylabel('kVar');
        subplot(2,3,5); plot(d.vvec(1:k,:,1)/1000);title('Charging rate: v, S=1, t=1'); ylabel('kW');
        subplot(2,3,6); plot(d.ugvec(1:k,:,1)/U0);title('U, G=1, t=1'); ylabel('Voltage, p.u.');
        
        figure(101), title('All times');
        subplot(2,3,1); plot(reshape(d.gluvec(1:k,1,:),k,nt));title('\lambda_U, S=1, t=1');
        subplot(2,3,2); plot(reshape(d.gllvec(1:k,1,:),k,nt));title('\lambda_L, S=1, t=1');
        subplot(2,3,4); plot(reshape(d.qvec(1:k,1,:),k,nt)/1000);title('Reactive support: q, S=1, t=1'); ylabel('kVar');
        subplot(2,3,5); plot(reshape(d.vvec(1:k,1,:),k,nt)/1000);title('Charging rate: v, S=1, t=1'); ylabel('kW');
        subplot(2,3,6); plot(reshape(d.ugvec(1:k,1,:),k,nt)/U0);title('U, G=1, t=1'); ylabel('Voltage, p.u.');
end

figure(105); plot(d.costtotal(1:k)); title('Total cost function');
figure(106); plot(sum(abs(d.subPower(1:k,:)),2)); title('Power consumption');
figure(107); plot(sum(abs(d.subLoss(1:k,:)),2)); title('Feeder''s total loss');
end

function saving(args,d,niter,k)
% global indent; 
if isempty(args.saveFn) || mod(k,args.plotFreq) > 0
    return;
end
if ~exist('tmp','dir'), mkdir('tmp'); end
save(['tmp/' args.saveFn '.mat'],'args.gs','d','niter');
end