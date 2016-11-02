function setVsource(o,field,val)
if isnumeric(val)
    val = num2str(val);
end
setval(o,'vsource.source',field,val);
end
