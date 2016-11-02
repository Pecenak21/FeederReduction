function wd = Convert2DSS_WireData(o)
% wd = Convert2DSS_WireData(o)
%
% PURPOSE : Converts wire data from object to OpenDSS
%
%
% INPUT :   o: struct that has all entries (e.g., o.R1)
%
% OUTPUT :  OpenDSS Wire Data object
%           see page 97 in OpenDSS manual V7.4.3(March 2012)

wd(length(o)) = dsswiredata;
for i = 1:length(o)
    wd(i).Name = o(i).ID;
    wd(i).Diam = o(i).Diameter;
    wd(i).GMRac = o(i).GMR;
    wd(i).Rdc = o(i).R25;
    wd(i).Rac = o(i).R50;
    wd(i).Normamps = o(i).Amps;
end


end

% 	obj.defaults = struct('Name','', ...
% 		'Rdc',[],... %dc Resistance, ohms per unit length (see Runits). Defaults to Rac if not specified.
% 		'Rac',[],... % Resistance at 60 Hz per unit length. Defaults to Rdc if not specified.
% 		'Runits','none',... % Length units for resistance: ohms per {mi|kft|km|m|Ft|in|cm } Default=none.
% 		'GMRac',[],...    % GMR at 60 Hz. Defaults to .7788*radius if not specified.
% 		'GMRunits','none',... % Units for GMR: {mi|kft|km|m|Ft|in|cm } Default=none.
% 		'Radius',[],... % Outside radius of conductor. Defaults to GMR/0.7788 if not specified.
% 		'Radunits','none',... % Units for outside radius: {mi|kft|km|m|Ft|in|cm } Default=none.
% 		'Normamps',[],... %Normal ampacity, amperes. Defaults to Emergency amps/1.5 if not specified.
% 		'Emergamps',[],... %Emergency ampacity, amperes. Defaults to 1.5 * Normal Amps if not specified.
% 		'Diam',[],... % Diam= Diameter; Alternative method for entering radius.
% 		'Like','');

% This class of data defines the raw conductor data that is used to compute the impedance for a
% line geometry.
% Note that you can use whatever units you want for any of the dimensional data – be sure to
% declare the units. Otherwise, the units are all assumed to match, which would be very rare for
% conductor data. Conductor data is usually supplied in a hodge?podge of units. Everything is
% converted to meters internally to the DSS
% Rdc= dc Resistance, ohms per unit length (see Runits). Defaults to Rac if not specified.
% Rac= Resistance at 60 Hz per unit length. Defaults to Rdc if not specified.
% Runits= Length units for resistance: ohms per {mi|kft|km|m|Ft|in|cm } Default=none.
% GMRac= GMR at 60 Hz. Defaults to .7788*radius if not specified.
% GMRunits= Units for GMR: {mi|kft|km|m|Ft|in|cm } Default=none.
% Radius= Outside radius of conductor. Defaults to GMR/0.7788 if not specified.
% Radunits= Units for outside radius: {mi|kft|km|m|Ft|in|cm } Default=none.
% Normamps= Normal ampacity, amperes. Defaults to Emergency amps/1.5 if not specified.
% Emergamps= Emergency ampacity, amperes. Defaults to 1.5 * Normal Amps if not specified.
% Diam= Diameter; Alternative method for entering radius.
% Like= Make like another object of this class: