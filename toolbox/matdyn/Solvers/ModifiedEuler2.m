function [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, t, eulerfailed, stepsize] = ModifiedEuler2(t, Xgen0, Pgen, Vgen0, Xexc0, Pexc, Vexc0, Xgov0, Pgov, Vgov0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, stepsize)

% [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, t,
% stepsize] = ModifiedEuler2(t, Xgen0, Pgen, Vgen0, Xexc0, Pexc, Vexc0,
% Xgov0, Pgov, Vgov0, Y, gbus, faulted_bus, theta0, genmodel, excmodel,
% govmodel, eulerfailed, stepsize) 
%
% Modified Euler ODE solver with check on interface errors
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Set up

eulerfailed = 0;

tol = 1e-8;
maxit = 20;

%% First Prediction Step
    
% EXCITERS
dFexc0 = Exciter(Xexc0, Pexc, Vexc0, excmodel);
Xexc_new = Xexc0 + stepsize.*dFexc0;

% GOVERNORS
dFgov0 = Governor(Xgov0, Pgov, Vgov0, govmodel);
Xgov_new = Xgov0 + stepsize.*dFgov0;
        
% GENERATORS
dFgen0 = Generator(Xgen0, Xexc_new, Xgov_new, Pgen, Vgen0, genmodel);
Xgen_new = Xgen0 + stepsize.*dFgen0;
       
Vexc_new = Vexc0;
Vgov_new = Vgov0;
Vgen_new = Vgen0;
    
for i=1:maxit
    Xexc_old = Xexc_new;
    Xgov_old = Xgov_new;
    Xgen_old = Xgen_new;
    
	Vexc_old = Vexc_new;
    Vgov_old = Vgov_new;
    Vgen_old = Vgen_new;
    
    % Calculate system voltages
    U_new = SolveNetwork(Xgen_new, Pgen, Ly, Uy, Py, gbus, genmodel);
    
    % Calculate machine currents and power
    [Id_new,Iq_new,Pe_new] = MachineCurrents(Xgen_new, Pgen, U_new(gbus), genmodel);
    
    % Update variables that have changed
    Vgen_new = [Id_new,Iq_new,Pe_new];
    Vexc_new = abs(U_new(gbus));
    Vgov_new = [Xgen_new(:,2)];    
    
    % Correct the prediction, and find new values of x
    % EXCITERS
    dFexc1 = Exciter(Xexc_old, Pexc, Vexc_new, excmodel);
    Xexc_new = Xexc0 + stepsize/2 .* (dFexc0 + dFexc1);
    
    % GOVERNORS
    dFgov1 = Governor(Xgov_old, Pgov, Vgov_new, govmodel);
    Xgov_new = Xgov0 + stepsize/2 .* (dFgov0 + dFgov1);
         
    % GENERATORS
    dFgen1 = Generator(Xgen_old, Xexc_new, Xgov_new, Pgen, Vgen_new, genmodel);
    Xgen_new = Xgen0 + stepsize/2 .* (dFgen0 + dFgen1);    
    
 
    
    % Calculate error
    Xexc_d = abs((Xexc_new-Xexc_old)');
    Xgov_d = abs((Xgov_new-Xgov_old)');
    Xgen_d = abs((Xgen_new-Xgen_old)');
    
	Vexc_d = abs((Vexc_new-Vexc_old)');
    Vgov_d = abs((Vgov_new-Vgov_old)');
    Vgen_d = abs((Vgen_new-Vgen_old)');
    
    errest = max( [max(max(Vexc_d)) max(max(Vgov_d)) max(max(Vgen_d)) max(max(Xexc_d)) max(max(Xgov_d)) max(max(Xgen_d)) ]);
       
    if errest < tol
        break    % solution found
    else
        if i==maxit
            U0 = U_new;
            Vexc0 = Vexc_new; Vgov0 = Vgov_new; Vgen0 = Vgen_new;
            Xgen0 = Xgen_new; Xexc0 = Xexc_new; Xgov0 = Xgov_new;
            Pgen0 = Pgen; Pexc0 = Pexc; Pgov0 = Pgov; 
            eulerfailed = 1;
            return;
        end
    end  
    
end


%% Update

U0 = U_new;
    
Vexc0 = Vexc_new;
Vgov0 = Vgov_new;
Vgen0 = Vgen_new;
    
Xgen0 = Xgen_new; 
Xexc0 = Xexc_new;
Xgov0 = Xgov_new;      
    
Pgen0 = Pgen; 
Pexc0 = Pexc;
Pgov0 = Pgov;    

return;