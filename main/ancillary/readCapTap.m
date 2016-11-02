function val = readCapTap(val)
% val is output of getval(o,circuit.transformer(1),'taps');
if iscell(val)
    val = cellfun(@readCapTapElmt, val,'uniformoutput',0);
else
    val = {readCapTapElmt(val)};
end
end

function val = readCapTapElmt(val)
val = strsplit(strtrim(val(2:end-1)),' ');
val = cell2mat(cellfun(@str2double,val,'uniformoutput',0));
end