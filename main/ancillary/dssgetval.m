function [ val ] = dssgetval( dssengine, element, valtype )
%DSSGETVAL(dssengine, element, valtype) 
%	return values of dssElement after simulation
%
% dssengine: dssEngine after simulation or first output of dssget. If the
% circuit struct is input instead, then automatically run dssget(dssengine) to get the
% simulation output.
%
% element: element(s) of interest either in dssObj or in string format: TYPE.NAME 
%         exp: c.fault(1) or 'fault.1'
%
% valtype: 
%			'currents' :(output depends on components)
%						1-phase load: 4 values of real and reactive currents of 1 phase with 2 terminals (in and out)
%						3-phase load: 8 values: first 6 are real and
%						reactive currents of 3 phases; last 2 are unknown
%			'voltages' in Magnitude and phase (degrees)
%			'powers'   in kw, kvar for all 3 phases
%			'dist'     in km
%           'ymatrix' or 'y' or 'ymat'  return ymatrix of the system
%
% examples:
%			o = dssget(c); v = dssgetval(o,c.load(1),'voltages')
%			v = dssgetval(c,c.generator,'dist')

if isstruct(dssengine)
	try
		o = dssget(dssengine);
	catch
		error('Incorrect input! Please check your input again!');
	end
else
	o = dssengine;
end

% handle dist separately 
if strcmp(valtype,'dist')
	% get the bus that the elem is connected to
	c = class(element);
	if ~strcmpi(c(1:3),'dss')
		error('Wrong element input for distance calculation! Input must be a dss object!');
	end
	switch c(4:end)
		case {'buslist'}
			error('Hasn''t implemented this feature for these components yet! Need to trace to the controlled components and find their locations!')
        case 'regcontrol'
            if ~isstruct(dssengine)
                error('Please input circuit (not simulation output) for this function to find distance of Voltage Regulators!');
            end
            [~,tr] = ismember(lower({element.transformer}),lower({dssengine.transformer.name}));
            for i = 1:length(element)
				b = dssengine.transformer(tr(i)).bus;
				bus{i} = b{1};
            end
        case 'capcontrol'
            if ~isstruct(dssengine)
                error('Please input circuit (not simulation output) for this function to find distance of Voltage Regulators!');
            end
            [~,tr] = ismember(lower({element.capacitor}),lower({dssengine.capacitor.name}));
            for i = 1:length(element)
				bus{i} = dssengine.capacitor(tr(i)).bus1;
            end
		case 'transformer'
			for i = 1:length(element)
				b = element(i).bus;
				bus{i} = b{1};
			end
        otherwise
            bus = {element.bus1};
	end
	bus = cleanBus(bus);
    
	% get dist of all buses
	bnm = o.ActiveCircuit.AllBusNames; %bus name
	bdist = o.ActiveCircuit.AllBusDistances'; % bus dist
	
	% get dist of desired elems based on the bus they connect to
	[a, b] = ismember(lower(bus),bnm);
	val = bdist(b(a));
	if any(a==0) 
		warning('Can''t find distance for some elements below:');
		disp(char(element(~a)));
    end
elseif strcmpi(valtype(1),'y') 
    ymat = o.activeCircuit.systemY;
    x = reshape(ymat,[2 length(ymat)/2]);
    si = sqrt(length(ymat)/2);
    x = x(1,:) - 1i*x(2,:);
    val = reshape(x,[si si])';
else
	if ~ischar(element) && strfind(class(element),'dss')
		for i = 1:length(element)
			c = class(element(i));
			el = [c(4:end) '.' element(i).Name];
			try
				val(i,:) = dssgetvalElem( el );
			catch
				vv = dssgetvalElem( el );
				a = numel(val(i-1,:));
				b = numel( vv );
				if a < b
					val(:,a+1:b) = nan;
				else
					vv(:,b+1:a) = nan;
				end
				val(i,:) = vv;
			end
		end
	else
		val = dssgetvalElem( element );
	end
end

function v = dssgetvalElem( elem)
	
	if strcmp(valtype, 'voltages')
		o.text.command = ['select ' elem];
		v = o.activeCircuit.ActiveElement.VoltagesMagAng;
	else
		o.text.command = ['select ' elem];
        % if property not found then try the 'getval' function 
        x = fixFieldname(valtype,o.activeCircuit.ActiveElement);
        if ~isempty(x)
            v = o.activeCircuit.ActiveElement.(x);
        else
            disp('Unknown Property!');
            v = [];
        end
	end
end

end

