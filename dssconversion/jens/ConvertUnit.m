function cf = ConvertUnit(Unit_impedance, Unit_length)
% cf = ConvertUnit(Unit_impedance, Unit_length)
%
% PURPOSE : gets conversion factor for converting ohms per length to ohms per unit length 
%
%
% INPUT :   Unit_impedance: OhmsPerMeter, OhmsPerCMeter, OhmsPerMile,
%                           OhmsPerFeet,OhmsPerCMeter
%           Unit_length: m, cm, mile, ft, kft
%
% OUTPUT : conversion factor (cf), multiply value given in 'Unit_impedance'
%          with cf to get Ohms per Unit_length,
%          e.g., 1 Ohms per mile, ft => cf = 0.000621371/3.28084
%                1 Ohms per mile * cf = 0.000621371/3.28084 Ohms per ft

if strcmp(Unit_impedance,'OhmsPerMeter')
    cf=1;
elseif strcmp(Unit_impedance,'OhmsPerKMeter')
    cf=1/1000;
elseif strcmp(Unit_impedance,'OhmsPerCMeter')
    cf=100;
elseif strcmp(Unit_impedance,'OhmsPerMile')
    cf=0.000621371;
elseif strcmp(Unit_impedance,'OhmsPerFeet')
    cf=3.28084;
elseif strcmp(Unit_impedance,'OhmsPerKFeet')
    cf=3.28084/1000;
else
    cf=1;
end
if strcmp(Unit_length,'m')
    cf=cf/1;
elseif strcmp(Unit_length,'cm')
    cf=cf/100;
elseif strcmp(Unit_length,'mile')
    cf=cf/0.000621371;
elseif strcmp(Unit_length,'ft')
    cf=cf/3.28084;
elseif strcmp(Unit_length,'kft')
    cf=cf/(3.28084/1000);
end


end