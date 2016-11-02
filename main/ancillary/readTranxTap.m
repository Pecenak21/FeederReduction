function val = readTranxTap(val)
% val is output of getval(o,circuit.transformer(1),'taps');
if iscell(val)
    val = cellfun(@readTranxTapElmt, val,'uniformoutput',0);
    val = cell2mat(val);
else
    val = readTranxTapElmt(val);
end
end

function val = readTranxTapElmt(val)
val = strsplit(val(2:end-1),', ');
val = cell2mat(cellfun(@str2double,val,'uniformoutput',0));
end