function c = convertPVModel(c,modelType)
% this function is to convert the pv system modelling style from complex
% model (needing Temp profile, etc.) to direct model using power output
% as input and vise versa.
%
% model 1: complex. Model 2: direct power input & output.

% default model to convert to is model 2 (direct power input & output)
if ~exist('modelType','var') || isempty(modelType)
    modelType = 2;
end

switch modelType
    case 2
        for i = 1:length(c.pvsystem)
            c.pvsystem(i).pmpp = c.pvsystem(i).kVA; % pv system output rating
            %c.pvsystem(i).kVA = c.pvsystem(i).kVA*; % inverter rating (make it 1.1*kw)
            c.pvsystem(i).irradiance = 1;
            % 3 phases
            %c.pvsystem(i).bus1 = cleanBus(c.pvsystem(i).bus1);
            %c.pvsystem(i).phases = 3;
            c.pvsystem(i).cutin = 0;
            c.pvsystem(i).cutout = 0;
            c.pvsystem(i).pf = [];
            c.pvsystem(i).kv = [];
        end
    otherwise
        error('This output model for PV system is not supported yet!')
end

end