function [y] = optGen(q,p,slim,T)
% finding optimal allowable generating state for generators
% q: reactive power
% p: active power
% slim: generation limits (0 as lower limit)
% T: number of time steps
% beta: max capacity of batw
cvx_begin quiet
   variable y(T)
   minimize( norm(q-y))
   subject to
		for t=1:T
			y(t) >= 0;
			y(t) <= sqrt(slim(t)^2 - p(t)^2); 
		end
cvx_end

if any(isnan(y))
    error;
end

end