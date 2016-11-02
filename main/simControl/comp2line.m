function c = comp2line(c,comp)
% replace circuit components by 1-cm default line
%
% input
%       c:        circuit
%       comp:     component
%
% example:
%       c = mvCompCloser(c,'transformer');

if strfind(class(comp),'dss')
    g = comp; type = strrep(class(comp),'dss','');
elseif ischar(comp)
    g = c.(comp); type = lower(g);
else
    error('not supported ''comp'' input!');
end

switch type
    case {'transformer','tranx'};
        for i = 1:length(g)
            l.Name = ['addedLine_' sprintf('%d',i)];
            l.bus1 = g(i).Buses{1};
            l.bus2 = g(i).Buses{2};
            l.length = 1;
            l.units = 'cm';
            c.line(end+1) = l;
        end
        if isfield(c,'regcontrol')
            rid = ismember(lower({c.regcontrol.transformer}),lower({g.Name}));
            c.regcontrol(rid) = [];
        end
        [~,gid] = ismember({g.Name},{c.transformer.Name});
        c.transformer(gid) = [];
    case {'vr','voltageregulator','regcontrol'}
        [~,tid] = ismember(lower({g.transformer}),lower({c.transformer.Name}));
        c = comp2line(c,c.transformer(tid));
    otherwise
        error('''comp'' type not supported. Please add code to support this comp');
end
end