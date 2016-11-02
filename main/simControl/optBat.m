function [y,z] = optBat(v,z0,vlim,t,beta)
% finding optimal allowable charging state 'y' for battery with
% v: array of charging states (desired charging rates)
% z0: initial capacity (1 for 100% full and 0 for empty bat)
% vlim: charging limits: 1 - lower limit, 2 - upper limit
% t: time
% beta: max capacity of bat (Wh)
nt = length(t);
dt = (t(2)-t(1))*24;% in hours
cvx_begin quiet
   variable y(nt) % bat charging rate (w)
   variable z(nt) % battery energy (wh)
   minimize( norm(v-y) ) 
   subject to
		z(1)==z0+y(1)/beta;
		z(1)<= 1;
		z(1)>= 0;
	    y(1) <= vlim(2);
		y(1) >= vlim(1);
		for t=2:nt
			z(t)==z(t-1)+y(t)/beta*dt;
			z(t)<= 1;
			z(t)>=0;
			y(t) <= vlim(2);
			y(t) >= vlim(1);
		end
cvx_end

if any(isnan(y))
    error('Optimization fails! Can''t find viable values for energy storage systems.');
end
end