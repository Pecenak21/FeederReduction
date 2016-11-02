function [Angles,Speeds,Eq_tr,Ed_tr,Efd,PM,Voltages,Stepsize,Errest,Time] = rundyn(casefile_pf, casefile_dyn, casefile_ev, mdopt)

% [Angles,Speeds,Eq_tr,Ed_tr,Efd,PM,Voltages,Stepsize,Errest,Time] =
% rundyn(casefile_pf, casefile_dyn, casefile_ev, mdopt)
% 
% Runs dynamic simulation
% 
% INPUTS
% casefile_pf = m-file with power flow data
% casefile_dyn = m-file with dynamic data
% casefile_ev = m-file with event data
% mdopt = options vector
% 
% OUTPUTS
% Angles = generator angles
% Speeds = generator speeds
% Eq_tr = q component of transient voltage behind reactance
% Ed_tr = d component of transient voltage behind reactance
% Efd = Excitation voltage
% PM = mechanical power
% Voltages = bus voltages
% Stepsize = step size integration method
% Errest = estimation of integration error
% Failed = failed steps
% Time = time points
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Begin timing
tic;
%% Add subdirectories to path

addpath([cd '/Solvers/']);
addpath([cd '/Models/Generators/']);
addpath([cd '/Models/Exciters/']);
addpath([cd '/Models/Governors/']);
addpath([cd '/Cases/Powerflow/']);
addpath([cd '/Cases/Dynamic/']);
addpath([cd '/Cases/Events/']);


%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% Options
if nargin < 4
	mdopt = Mdoption;   
end
method = mdopt(1);
tol = mdopt(2);
minstepsize = mdopt(3);
maxstepsize = mdopt(4);
output = mdopt(5);
plots = mdopt(6);


%% Load all data

% Load dynamic simulation data
if output; disp('> Loading dynamic simulation data...'); end
global freq
[freq,stepsize,stoptime] = Loaddyn(casefile_dyn);

% Load generator data
Pgen0 = Loadgen(casefile_dyn, output);

% Load exciter data
Pexc0 = Loadexc(casefile_dyn);

% Load governor data
Pgov0 = Loadgov(casefile_dyn);

% Load event data
if ~isempty(casefile_ev)
    [event,buschange,linechange] = Loadevents(casefile_ev);
else
    event=[];
end

genmodel = Pgen0(:,1);
excmodel = Pgen0(:,2);
govmodel = Pgen0(:,3);


%% Initialization: Power Flow 

% Power flow options
mpopt=mpoption;
% mpopt(31)=0;
% mpopt(32)=0;

% Run power flow
[baseMVA, bus, gen, branch, success] = runpf(casefile_pf,mpopt);
if ~success
    fprintf('> Error: Power flow did not converge. Exiting...\n')
    return;
else
    if output; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> Power flow converged\n'); end
end

U0=bus(:,VM).*(cos(bus(:,VA)*pi/180) + j.*sin(bus(:,VA)*pi/180));
U00=U0;
% Get generator info
on = find(gen(:, GEN_STATUS) > 0);     %% which generators are on?
gbus = gen(on, GEN_BUS);               %% what buses are they at?
ngen = length(gbus);

nbus = length(U0);


%% Construct augmented Ybus 
if output; disp('> Constructing augmented admittance matrix...'); end
Pl=bus(:,PD)./baseMVA;                  %% load power
Ql=bus(:,QD)./baseMVA;

xd_tr = zeros(ngen,1);
xd_tr(genmodel==2) = Pgen0(genmodel==2,8); % 4th order model: xd_tr column 8
xd_tr(genmodel==1) = Pgen0(genmodel==1,7); % classical model: xd_tr column 7

[Ly, Uy, Py] = AugYbus(baseMVA, bus, branch, xd_tr, gbus, Pl, Ql, U0);


%% Calculate Initial machine state
if output; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> Calculating initial state...\n'); end
[Efd0, Xgen0] = GeneratorInit(Pgen0, U0(gbus), gen, baseMVA, genmodel);

omega0 = Xgen0(:,2);

[Id0,Iq0,Pe0] = MachineCurrents(Xgen0, Pgen0, U0(gbus), genmodel);
Vgen0 = [Id0, Iq0, Pe0];


%% Exciter initial conditions

Vexc0 = [abs(U0(gbus))];
[Xexc0,Pexc0] = ExciterInit(Efd0, Pexc0, Vexc0, excmodel);


%% Governor initial conditions

Pm0 = Pe0;
[Xgov0, Pgov0] = GovernorInit(Pm0, Pgov0, omega0, govmodel);
Vgov0 = omega0;


%% Check Steady-state

Fexc0 = Exciter(Xexc0, Pexc0, Vexc0, excmodel);
Fgov0 = Governor(Xgov0, Pgov0, Vgov0, govmodel);
Fgen0 = Generator(Xgen0, Xexc0, Xgov0, Pgen0, Vgen0, genmodel);

% Check Generator Steady-state
if sum(sum(abs(Fgen0))) > 1e-6
    fprintf('> Error: Generator not in steady-state\n> Exiting...\n')
    return;
end
% Check Exciter Steady-state
if sum(sum(abs(Fexc0))) > 1e-6
	fprintf('> Error: Exciter not in steady-state\n> Exiting...\n')
	return;
end
% Check Governor Steady-state
if sum(sum(abs(Fgov0))) > 1e-6
    fprintf('> Error: Governor not in steady-state\n> Exiting...\n')
    return;
end

if output; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> System in steady-state\n'); end


%% Initialization of main stability loop
t=-0.02; % simulate 0.02s without applying events
errest=0;
failed=0;
eulerfailed = 0;

if method==3 || method==4
    stepsize = minstepsize;
end


if ~output
    fprintf('                   ')
end

ev=1;
eventhappened = false;
i=0;


%% Allocate memory for variables

if output; fprintf('> Allocate memory..'); end
chunk = 5000;

Time = zeros(chunk,1); Time(1,:) = t;
Errest = zeros(chunk,1); Errest(1,:) = errest;
Stepsize = zeros(chunk,1); Stepsize(1,:) = stepsize;

% System variables
Voltages = zeros(chunk, length(U0)); Voltages(1,:) = U0.';

% Generator
Angles = zeros(chunk,ngen); Angles(1,:) = Xgen0(:,1).*180./pi;
Speeds = zeros(chunk,ngen); Speeds(1,:) = Xgen0(:,2)./(2.*pi.*freq);
Eq_tr = zeros(chunk,ngen); Eq_tr(1,:) = Xgen0(:,3);
Ed_tr = zeros(chunk,ngen); Ed_tr(1,:) = Xgen0(:,4);

% Exciter and governor
Efd = zeros(chunk,ngen); Efd(1,:) = Efd0(:,1);  
PM = zeros(chunk,ngen); PM(1,:) = Pm0(:,1);


%% Main stability loop
while t < stoptime + stepsize

    %% Output    
    i=i+1;
    if mod(i,45)==0 && output
       fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> %6.2f%% completed', t/stoptime*100)
    end  
    
    %% Numerical Method
    switch method
        case 1
            [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, t, newstepsize] = ModifiedEuler(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, stepsize);            
        case 2
            [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, t, newstepsize] = RungeKutta(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, stepsize);
        case 3
            [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, errest, failed, t, newstepsize] = RungeKuttaFehlberg(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, tol, maxstepsize, stepsize);
        case 4          
            [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, errest, failed, t, newstepsize] = RungeKuttaHighamHall(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, tol, maxstepsize, stepsize);                 
        case 5
            [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, t, eulerfailed, newstepsize] = ModifiedEuler2(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, stepsize);
    end
        
    if eulerfailed
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> Error: No solution found. Try lowering tolerance or increasing maximum number of iterations in ModifiedEuler2. Exiting... \n')
    	return;
    end        
    
	if failed
        t = t-stepsize;
    end

    % End exactly at stop time
    if t + newstepsize > stoptime
    	newstepsize = stoptime - t;
    elseif stepsize < minstepsize
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> Error: No solution found with minimum step size. Exiting... \n')
    	return;
    end;
    
    
    %% Allocate new memory chunk if matrices are full
    if i>size(Time,1)
        Stepsize = [Stepsize; zeros(chunk,1)];Errest = [Errest; zeros(chunk,1)];Time = [Time; zeros(chunk,1)];
        Voltages = [Voltages; zeros(chunk,length(U0))];Efd = [Efd; zeros(chunk,ngen)];PM = [PM; zeros(chunk,ngen)];
        Angles=[Angles;zeros(chunk,ngen)];Speeds=[Speeds;zeros(chunk,ngen)];Eq_tr=[Eq_tr;zeros(chunk,ngen)];Ed_tr=[Ed_tr;zeros(chunk,ngen)];
    end

    
    %% Save values
    Stepsize(i,:) = stepsize.';
    Errest(i,:) = errest.';
    Time(i,:) = t;
    
    Voltages(i,:) = U0.';

    % exc
    Efd(i,:) = Xexc0(:,1).*(genmodel>1); % Set Efd to zero when using classical generator model  
    
    % gov
    PM(i,:) = Xgov0(:,1);
   
    % gen
	Angles(i,:) = Xgen0(:,1).*180./pi;
    Speeds(i,:) = Xgen0(:,2)./(2.*pi.*freq);
    Eq_tr(i,:) = Xgen0(:,3);
    Ed_tr(i,:) = Xgen0(:,4);   
    
    
    %% Adapt step size if event will occur in next step
    if ~isempty(event) && ev <= size(event,1) && (method == 3 || method == 4)      
        if t + newstepsize >= event(ev,1)
            if event(ev,1) - t < newstepsize
                newstepsize = event(ev,1) - t;
            end
        end
    end 

    
    %% Check for events
    if ~isempty(event) && ev <= size(event,1)    
           
        for k=ev:size(event,1) % cycle through all events ..   
            if abs(t-event(ev,1))>10*eps ||  ev > size(event,1) %.. that happen on time t               
                break;
            else
                eventhappened = true;
            end

                switch event(ev,2)
                    case 1
                        bus(buschange(ev,2),buschange(ev,3)) = buschange(ev,4);

                    case 2
                        branch(linechange(ev,2),linechange(ev,3)) = linechange(ev,4);                        
                end 
                ev=ev+1;
        end          
            
        if eventhappened
            % Refactorise
            [Ly, Uy, Py] = AugYbus(baseMVA, bus, branch, xd_tr, gbus, bus(:,PD)./baseMVA, bus(:,QD)./baseMVA, U00);
            U0 = SolveNetwork(Xgen0, Pgen0, Ly, Uy, Py, gbus, genmodel);            
            
            [Id0,Iq0,Pe0] = MachineCurrents(Xgen0, Pgen0, U0(gbus), genmodel);
            Vgen0 = [Id0,Iq0,Pe0];
            Vexc0 = abs(U0(gbus));
            
            % decrease stepsize after event occured
            if method==3 || method==4
                newstepsize = minstepsize;
            end
            
            i=i+1; % if event occurs, save values at t- and t+

            %% Save values
            Stepsize(i,:) = stepsize.';
            Errest(i,:) = errest.';            
            Time(i,:) = t;
    
            Voltages(i,:) = U0.';

            % exc
            Efd(i,:) = Xexc0(:,1).*(genmodel>1); % Set Efd to zero when using classical generator model  
    
            % gov
            PM(i,:) = Xgov0(:,1);
   
            % gen
            Angles(i,:) = Xgen0(:,1).*180./pi;
            Speeds(i,:) = Xgen0(:,2)./(2.*pi.*freq);
            Eq_tr(i,:) = Xgen0(:,3);
            Ed_tr(i,:) = Xgen0(:,4);
            
            eventhappened = false;
        end
    end    

    
    %% Advance time    
    stepsize = newstepsize;
    t = t + stepsize;
    
end % end of main stability loop


%% Output
if output
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> 100%% completed')
else
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
end
simulationtime=toc;
if output; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> Simulation completed in %5.2f seconds\n', simulationtime); end

% Save only the first i elements
Angles = Angles(1:i,:);
Speeds = Speeds(1:i,:);
Eq_tr = Eq_tr(1:i,:);
Ed_tr = Ed_tr(1:i,:);

Efd = Efd(1:i,:);
PM = PM(1:i,:);

Voltages = Voltages(1:i,:);

Stepsize = Stepsize(1:i,:);
Errest = Errest(1:i,:);
Time = Time(1:i,:);

%% Clean up
rmpath([cd '/Solvers/']);
rmpath([cd '/Models/Generators/']);
rmpath([cd '/Models/Exciters/']);
rmpath([cd '/Models/Governors/']);
rmpath([cd '/Cases/Powerflow/']);
rmpath([cd '/Cases/Dynamic/']);
rmpath([cd '/Cases/Events/']);

%% Plot
if plots
    
close all
    
figure
xlabel('Time [s]')
ylabel('Angle [deg]')
hold on
plot(Time,Angles)
axis([0 Time(end) -1 1])
axis 'auto y'

figure
xlabel('Time [s]')
ylabel('Speed [pu]')
hold on
plot(Time,Speeds)
axis([0 Time(end) -1 1])
axis 'auto y'

figure
xlabel('Time [s]')
ylabel('Voltage [pu]')
hold on
plot(Time,abs(Voltages))
axis([0 Time(end) -1 1])
axis 'auto y'

figure
xlabel('Time [s]')
ylabel('Excitation voltage [pu]')
hold on
plot(Time,Efd)
axis([0 Time(end) -1 1])
axis 'auto y'

figure
xlabel('Time [s]')
ylabel('Turbine Power [pu]')
hold on
plot(Time,PM)
axis([0 Time(end) -1 1])
axis 'auto y'

figure
hold on
xlabel('Time [s]')
ylabel('Step size')
plot(Time,Stepsize,'-o')
axis([0 Time(end) -1 1])
axis 'auto y'

end

return;