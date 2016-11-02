function [ res ] = getLoadsDownstream( c,firstbus )
%GETBUSDOWNSTREAM gives you: a list of  1. downstream buses (cleaned)
%                                       2. downstream bus numbers
%                                       3. downstream line numbers 
%                                       4. downstream load numbers
%                                       5. downstream load buses (with phases)
%from the number of the bus from which you want to start and the circuit.
%The function also plot the feeder to let you check (downstream load 
%objects are in the field "VviolationMax").
%
%   Inputs:
%       - "c":              circuit object
%       - "firstbus":       number (in buslist object) of the bus from 
%                                   which the function will start
%   Output:
%       - "res":            structure with all results in fields:
%                                       1. buses
%                                       2. numbuslist
%                                       3. linenum
%                                       4. loadnum
%                                       5. loadbus

% get buslist
res=struct('buses',{{''}},'numbuslist',[1],'linenum',[1],'loadnum',[],'loadbus',{{''}});
counter=1;
while counter<length(res.buses)+1
    aa = find(ismember(cleanBus({c.line.bus1}),c.buslist.id(firstbus)));
    for i=1:length(aa)
        res.linenum(end+1) = aa(i);
        res.numbuslist(end+1) = find(ismember(c.buslist.id,cleanBus(c.line(aa(i)).bus2)));
        res.buses(end+1) = c.buslist.id(res.numbuslist(end));
    end
    counter=counter+1;
    if (counter<length(res.buses)+1)
        firstbus=res.numbuslist(counter);
    end
end
res.linenum=unique(res.linenum);
res.numbuslist = unique(res.numbuslist);
res.buses=unique(res.buses);

% get loadlist
for j =1:length(res.buses)-1
    if sum(ismember(cleanBus({c.load.bus1}),res.buses{j+1}))>0
        aa=find(ismember(cleanBus({c.load.bus1}),res.buses{j+1}));
        for l=1:length(aa)
            res.loadnum(end+1) = aa(l);
            res.loadbus{end+1} = c.load(aa(l)).bus1;
        end
    end
end
    

% plot it to check
for j =1:length(res.loadbus)-1
    c.VviolationMax(j) =dsscapacitor;
    c.VviolationMax(j).Name = res.loadbus{j+1};
    c.VviolationMax(j).Bus1 = res.loadbus{j+1};
end
h = circuitVisualizer(c);
set(h,'Name',[get(h,'Name') ' downstream loads shown as ''VviolationMax'' from specified bus']);
end

