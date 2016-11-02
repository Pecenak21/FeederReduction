function c = rmCapacitor(c)
if isfield(c,'capcontrol')
    c = rmfield(c,'capcontrol');
end
if isfield(c,'capacitor')
    c = rmfield(c,'capacitor');
end
end