%% test circuit with 2 nodes, 1 of them is the PCC
% circuit (sourcebus)
x = dsscircuit;
x.Name = 'lacopf';
x.bus1 = 'b0';
x.Phases = 1;
c.circuit = x;

% line
R = .2; X = 1; L = 1;
x = dssline;
x.Name = 'line1';
x.Phases = 1;
x.length = L;
x.rmatrix = [R];
x.xmatrix = [X];
x.bus1 = 'b0';
x.bus2 = 'b1';
c.line = x;

% load
x = dssload;
x.Name = 'load1';
x.model = 1;
x.kW = 100;
x.kvar = 10;
x.bus1 = 'b1';
x.phases = 1;
c.load = x;

% run simulation
o = dssget(c);
ymat = dssgetval(o,'','y');
y = calculateYMatrix({'b0','b1'},c);